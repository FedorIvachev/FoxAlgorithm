#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

int ProcRank = 0;     // Rank of current process
int GridSize;         // Size of virtual processor grid
int GridCoords[2];    // Coordinates of current processor in grid
MPI_Comm GridComm;    // Grid communicator
MPI_Comm ColComm;     // Column communicator
MPI_Comm RowComm;     // Row communicator



void MatrixInit1(int *A, int n){
    srand(time(NULL));
    int num = rand();
    int i, j;
    for (i = 0; i < n; ++i){
        for (j = 0; j < n; ++j){
            A[i * n + j] = num;
            num = rand();
        }
    }
}


void MatrixInit(int *A, int n){
    int num = 1;
    int i, j;
    for (i = 0; i < n; ++i){
        for (j = 0; j < n; ++j){
            A[i * n + j] = num;
            num++;
        }
    }
}

void PrintMatrix(int *A, int n){
    int i, j;
    for (i = 0; i < n; ++i){
        for (j = 0; j < n; ++j){
            printf("%d ", A[i * n + j]);
        }
        printf("\n");
    }
}


void CreateGridCommunicators() {
    int DimSize[2];  // Number of processes in each dimension of the grid
    int Periodic[2]; // =1, if the grid dimension should be periodic
    int Subdims[2];  // =1, if the grid dimension should be fixed
    DimSize[0] = GridSize;
    DimSize[1] = GridSize;
    Periodic[0] = 0;
    Periodic[1] = 0;    // Creation of the Cartesian communicator
    MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize, Periodic, 1, &GridComm);
    // Determination of the cartesian coordinates for every process
    MPI_Cart_coords(GridComm, ProcRank, 2, GridCoords);
    // Creating communicators for rows
    Subdims[0] = 0;  // Dimensionality fixing
    Subdims[1] = 1;  // The presence of the given dimension in the subgrid
    MPI_Cart_sub(GridComm, Subdims, &RowComm);
    // Creating communicators for columns
    Subdims[0] = 1;
    Subdims[1] = 0;
    MPI_Cart_sub(GridComm, Subdims, &ColComm);

}



void Distribute(int *Matrix, int *Block, int Size, int BlockSize)
{
    int *Buff = new int[BlockSize * Size];// буфер для рассылки блоков
    int i;
    if (GridCoords[1] == 0) 
        MPI_Scatter(Matrix, BlockSize*Size, MPI_INT, Buff, BlockSize * Size, MPI_INT, 0, ColComm);
        //распределение строк и столбцов между первыми проц-ми 
    for (i = 0; i < Size; ++i) {
        MPI_Scatter(&Buff[i * Size], BlockSize, MPI_INT, &Block[i * BlockSize], BlockSize, MPI_INT, 0, RowComm);
        //дальнейшее распределение строк между процессами
    }
    delete [] Buff;
}

void DataDistribution(int *pAMatrix, int *pBMatrix, int *pMatrixAblock, int *pBblock, int Size, int BlockSize){
    Distribute(pAMatrix, pMatrixAblock, Size, BlockSize);
    Distribute(pBMatrix, pBblock, Size, BlockSize);
}

void ResultCollection(int *pCMatrix, int *pCblock, int Size, int BlockSize) {
    int *Buff = new int[BlockSize * Size];// буфер для рассылки блоков
    int i;
    for (i = 0; i < Size; ++i){   // сборка матрицы С из блоков
        MPI_Gather(&pCblock[i * BlockSize], BlockSize, MPI_INT, &Buff[i * Size], BlockSize, MPI_INT, 0, RowComm);
    }

    if (GridCoords[1] == 0) {   // сборка из строк  
        MPI_Gather(Buff, BlockSize * Size, MPI_INT, pCMatrix, BlockSize * Size, MPI_INT, 0, ColComm);
    } 

    delete [] Buff;
} 


void ABlockCommunication (int iter, int *pAblock, int* pMatrixAblock, int BlockSize) {
    // Defining the leading process of the process grid row
    int Pivot = (GridCoords[0] + iter) % GridSize;
    // Copying the transmitted block in a separate memory buffer
    if (GridCoords[1] == Pivot) {
        for (int i = 0; i < BlockSize*BlockSize; ++i)
            pAblock[i] = pMatrixAblock[i];
    }
    // Block broadcasting
    MPI_Bcast(pAblock, BlockSize*BlockSize, MPI_INT, Pivot, RowComm);
}


void BblockCommunication (int *pBblock, int BlockSize) {
    MPI_Status Status;
    int NextProc = GridCoords[0] + 1;
    if ( GridCoords[0] == GridSize-1 ) NextProc = 0;
    int PrevProc = GridCoords[0] - 1;
    if ( GridCoords[0] == 0 ) PrevProc = GridSize-1;
    MPI_Sendrecv_replace( pBblock, BlockSize*BlockSize, MPI_INT,
                         NextProc, 0, PrevProc, 0, ColComm, &Status);
}


void BlockMultiplication (int *pAblock, int *pBblock, int *pCblock, int BlockSize) {   
// вычисление произведения матричных блоков   
    int i, j, k, temp;
    for (i = 0; i < BlockSize; i++) {     
        for (j = 0; j < BlockSize; j++) {       
            temp = 0;       
            for (k = 0; k < BlockSize; k++ )         
                temp += pAblock [i*BlockSize + k] * pBblock [k*BlockSize + j];    
            pCblock [i * BlockSize + j] += temp;     
        }   
    } 
}

void ParallelResultCalculation(int* pAblock, int* pMatrixAblock,
                               int* pBblock, int* pCblock, int BlockSize) {
    for (int iter = 0; iter < GridSize; iter ++) {
        // Sending blocks of matrix A to the process grid rows
        ABlockCommunication (iter, pAblock, pMatrixAblock, BlockSize);
        // Block multiplication
        BlockMultiplication(pAblock, pBblock, pCblock, BlockSize);
        // Cyclic shift of blocks of matrix B in process grid columns
        BblockCommunication(pBblock, BlockSize);
    }
}

int main (int argc, char * argv[] ) {
    int* pAMatrix;  // The first argument of matrix multiplication
    int* pBMatrix;  // The second argument of matrix multiplication
    int* pCMatrix;  // The result matrix
    int Size;       // Size of matricies
    int BlockSize;  // Sizes of matrix blocks on current process
    int *pAblock;   // Initial block of matrix A on current process
    int *pBblock;   // Initial block of matrix B on current process
    int *pCblock;   // Block of result matrix C on current process
    int ProcNum;
    int i;
    double timerMpi;
    int *pMatrixAblock;
    int Start, Finish, Duration;
    //setvbuf(stdout, 0, _IONBF, 0);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    GridSize = sqrt((int)ProcNum);
    if (ProcNum != GridSize*GridSize) {
        if (ProcRank == 0) {
            printf ("Number of processes must be a perfect square \n");
        }
    } else {
        if (ProcRank == 0)
            printf("Parallel matrix multiplication program\n");
        // Creating the grid, row and column communcators
        CreateGridCommunicators();
        // Memory allocation and initialization of matrix elements
        int i;
        if (ProcRank == 0) {
            do {
                printf("\nEnter size of the initial objects: ");
                scanf("%d", &Size);
                if (Size%GridSize != 0) {
                    printf ("Size of matricies must be divisible by the grid size! \n");
                }
            }
            while (Size%GridSize != 0);
        }
        MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);

        BlockSize = Size/GridSize;
        pAblock = new int [BlockSize * BlockSize];
        pBblock = new int [BlockSize * BlockSize];
        pCblock = new int [BlockSize * BlockSize];
        pMatrixAblock = new int [BlockSize * BlockSize * sizeof(int)];
        for (i = 0; i < BlockSize*BlockSize; ++i) {
            pCblock[i] = 0;
        }
        if (ProcRank == 0) {
            pAMatrix = new int [Size * Size];
            pBMatrix = new int [Size * Size];
            pCMatrix = new int [Size * Size];
            MatrixInit(pAMatrix, Size);
            MatrixInit(pBMatrix, Size);
            PrintMatrix(pAMatrix, Size);
            PrintMatrix(pBMatrix, Size);
        }
        DataDistribution(pAMatrix, pBMatrix, pMatrixAblock, pBblock, Size,
                         BlockSize);
        // Execution of Fox method
        ParallelResultCalculation(pAblock, pMatrixAblock, pBblock, pCblock, BlockSize);
        ResultCollection(pCMatrix, pCblock, Size, BlockSize);
        //MPI_Reduce(pCMatrix, pCMatrix, Size * Size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(GridComm);

        PrintMatrix(pCMatrix, Size);
        delete [] pAblock;
        delete [] pBblock;
        delete [] pCblock;
        delete [] pMatrixAblock;


        if (ProcRank == 0){
            delete [] pAMatrix;
            delete [] pBMatrix;
            delete [] pCMatrix;
        }

        
        // Нужно дописать завершение процессов
    }
    printf("%d\n", ProcRank);

    MPI_Finalize();
    return 0;
}


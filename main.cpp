#include <iostream>
#include <stdlib.h>
#include <mpi.h>
#include <ctime>
#include <math.h>

#define SZ 16

int MatrixA[SZ][SZ];
int MatrixB[SZ][SZ];
int MatrixC[SZ*SZ];

int ProcRank;
int GridSize;
int GridCoords[2];
MPI_Comm GridComm;
MPI_Comm ColComm;
MPI_Comm RowComm;


void CreateGridComm();

void BlockMult(int **a, int **b, int **c, int size);

void FoxCalculation(int Size, int **a, int **b, int **c);

void VecToMatrix(int *Buff, int **a, int size);

void MatrixToVec(int *Buff, int **a, int size);

int main(int argc, char *argv[]) {
    int i, j;
    int BlockSize;
    int ProcNum;
    int **CurA;
    int **CurB; 
    int **CurC;
    double b, e;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    GridSize = sqrt((int)ProcNum);
    BlockSize = SZ / GridSize;

    CreateGridComm();
    for (i = 0; i < SZ; ++i) {
        for (j = 0; j < SZ; ++j) {
            MatrixA[i][j] = rand() % 70;
            MatrixB[i][j] = rand() % 70;
        }
    }

    int base_row = GridCoords[0] * BlockSize;
    int base_col = GridCoords[1] * BlockSize;

    CurA = new int*[SZ];
    CurB = new int*[SZ];
    CurC = new int*[SZ];

    for(i = 0; i < SZ; ++i)
    {
        CurA[i] = new int[SZ];
        CurB[i] = new int[SZ];
        CurC[i] = new int[SZ];
    }

    for (i = base_row; i < base_row + BlockSize; ++i) {
        for (j = base_col; j < base_col + BlockSize; ++j) {
            CurA[i - (base_row)][j - (base_col)] = MatrixA[i][j];
            CurB[i - (base_row)][j - (base_col)] = MatrixB[i][j];
            CurC[i - (base_row)][j - (base_col)] = 0;
        }
    }

    MPI_Barrier(GridComm);

    if (ProcRank == 0){
        b = MPI_Wtime();
    }

    FoxCalculation(SZ, CurA, CurB, CurC);

    MPI_Barrier(GridComm);

    if (ProcRank == 0)
    {
        e = MPI_Wtime();
        std::cout << e-b << std::endl;
    }

    int *VecC = new int[SZ * SZ];
    int *CurBuff = new int[BlockSize * BlockSize];
    MatrixToVec(CurBuff, CurC, BlockSize);

    MPI_Gather(CurBuff, BlockSize * BlockSize, MPI_INT, VecC, BlockSize * BlockSize, MPI_INT, 0, GridComm);
    MPI_Barrier(GridComm);

    if (ProcRank == 0) {
        int k = 0;
        int i1, j1;
        for (i1 = 0; i1 < GridSize; ++i1) {
            for (j1 = 0; j1 < GridSize; ++j1) {
                for (i = i1 * BlockSize; i < i1 * BlockSize + BlockSize; ++i) {
                    for (j = j1 * BlockSize; j < j1 * BlockSize + BlockSize; ++j) {
                        MatrixC[i * SZ + j] = VecC[k];
                        k++;
                    }
                }
            }
        }
    }
    MPI_Finalize();
    return 0;
}   

void CreateGridComm() {
    int Dims[2];
    int Per[2];
    int SubDims[2];
    Dims[0] = GridSize;
    Dims[1] = GridSize;
    Per[0] = 0;
    Per[1] = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, Dims, Per, 1, &(GridComm));
    MPI_Cart_coords(GridComm, ProcRank, 2, GridCoords);
    SubDims[0] = 0;
    SubDims[1] = 1;
    MPI_Cart_sub(GridComm, SubDims, &(RowComm));
    SubDims[0] = 1;
    SubDims[1] = 0;
    MPI_Cart_sub(GridComm, SubDims, &(ColComm));
}

void FoxCalculation(int Size, int **a, int **b, int **c) {
    int **TempBlockA;
    int *Buff;
    int i, stage;
    int CurRoot; // Leading process
    int BlockSize, PrevProcess, NextProcess;
    MPI_Status status;
    BlockSize = Size / GridSize;

    PrevProcess = (GridCoords[0] + 1) % GridSize;
    NextProcess = (GridCoords[0] + GridSize - 1) % GridSize;

    // init of Buff
    Buff = new int[BlockSize * BlockSize];
    for (int i = 0; i < BlockSize * BlockSize; ++i)
        Buff[i] = 0;


    // init of TempBlockA
    TempBlockA = new int*[BlockSize];
    for(int i = 0; i < BlockSize; ++i)
        TempBlockA[i] = new int[BlockSize];
    for (int i = 0; i < BlockSize; ++i)
        for (int j = 0; j < BlockSize; ++j)
            TempBlockA[i][j] = 0;


    for (stage = 0; stage < GridSize; stage++) {
        CurRoot = (GridCoords[0] + stage) % GridSize;
        if (GridCoords[1] == CurRoot) {
            MatrixToVec(Buff, a, BlockSize);
            MPI_Bcast(Buff, BlockSize * BlockSize, MPI_INT, CurRoot, RowComm);
            VecToMatrix(Buff, a, BlockSize);
            BlockMult(a, b, c, BlockSize);
        } else {
            MatrixToVec(Buff, TempBlockA, BlockSize);
            MPI_Bcast(Buff, BlockSize * BlockSize, MPI_INT, CurRoot, RowComm);
            VecToMatrix(Buff, TempBlockA, BlockSize);
            BlockMult(TempBlockA, b, c, BlockSize);
        }
        MatrixToVec(Buff, b, BlockSize);
        MPI_Sendrecv_replace(Buff, BlockSize * BlockSize, MPI_INT, NextProcess, 0, PrevProcess, 0, ColComm, &status);
        VecToMatrix(Buff, b, BlockSize);
    }
}

void BlockMult(int **A, int **B, int **C, int size) {
    int i, j, k;
    int sum = 0;
    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) {
            sum = 0;
            for (k = 0; k < size; k++)
                sum += (A[i][k] * B[k][j]);
            C[i][j] += sum;
        }
    }
}

void VecToMatrix(int *vec, int **matrix, int size) {
    int i, j;
    int cur = 0;
    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) {
            matrix[i][j] = vec[cur];
            cur++;
        }
    }
}

void MatrixToVec(int *vec, int **matrix, int size) {
    int i, j;
    int cur = 0;
    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) {
            vec[cur] = matrix[i][j];
            cur++;
        }
    }
}

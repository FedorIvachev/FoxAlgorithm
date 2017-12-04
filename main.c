#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>


void PrintMatrix(int *a, int n){
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%d ", a[i*n+j]);
        }
        printf("\n");
    }
}

int CheckAnswer(int *a, int *b, int n)
{
    int i, j, k;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (a[i * n + j] != b[i * n + j]){
                return 0;
            }
        }
    }
    return 1;
}

void BasicMultiplication (int *A, int *B, int *C, int n)
{   
    int i, j, k;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            int sum = 0;
            for (k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] += sum;
        }
    }
}

void FoxAlgorythm (int *A, int *B, int *C, int n, int nProcNum){
    int stage;
    #pragma omp parallel private(stage) shared(A, B, C) num_threads(nProcNum)
    {
        int i, j, k, i1, j1;
        int GridSize = sqrt(nProcNum);
        //printf("GridSize: %d\n", GridSize);
        int PrNum = omp_get_thread_num();
        //printf("PrNum: %d\n", PrNum);
        i1 = PrNum / GridSize;
        j1 = PrNum % GridSize;
        int *A1;
        int *B1;
        int *C1;
        for (stage = 0; stage < GridSize; stage++){
            A1 = A + (n * i1 + ((i1 + stage) % GridSize)) * (n/GridSize);
            B1 = B + (n * ((i1 + stage) % GridSize) + j1) * (n/GridSize);
            C1 = C + (n * i1 + j1) * (n/GridSize);
            for (i = 0; i < n / GridSize; i++){
                for (j = 0; j < n / GridSize; j++){
                    for (k = 0; k < n / GridSize; k++){
                        C1[i*n+j] += A1[i*n+k] * B1[k*n+j];                    
                    }
                }                           
            }
        }
    }
}



int main(int argc, char *argv[])
{
    int mSize, ProcNum;

    sscanf(argv[1], "%d", &mSize);
    sscanf(argv[2], "%d", &ProcNum);
    
    int *pAMatrix = malloc(mSize * mSize * sizeof(int));
    int *pBMatrix = malloc(mSize * mSize * sizeof(int));
    int *woOpenMp = malloc(mSize * mSize * sizeof(int));
    int *wOpenMp = malloc(mSize * mSize * sizeof(int));
    
    int num = 1;
    int i, j;
    for (i = 0; i < mSize; i++) {
        for (j = 0; j < mSize; j++) {
            pAMatrix[i * mSize + j] = num;
            pBMatrix[i * mSize + j] = num;
            num++;
            woOpenMp[i * mSize + j] = 0;
            wOpenMp[i * mSize + j] = 0;
        }
    }
    double timerOpenMp = omp_get_wtime();
    FoxAlgorythm(pAMatrix, pBMatrix, wOpenMp, mSize, ProcNum);
    //BasicMultiplication(pAMatrix, pBMatrix, wOpenMp, mSize); printf("Basic alg: \n");

    timerOpenMp = omp_get_wtime() - timerOpenMp;

    printf("%.10lf ||||||| ProcNum: %d mSize: %d \n", timerOpenMp, ProcNum, mSize);
    
    //check answer
    /*
    PrintMatrix(wOpenMp, mSize);
    BasicMultiplication(pAMatrix, pBMatrix, woOpenMp, mSize);
    PrintMatrix(woOpenMp, mSize);
    if (CheckAnswer(woOpenMp, wOpenMp, mSize) == 1){
        printf("Answer is correct\n");
    }*/

    free(pAMatrix);
    free(pBMatrix);
    free(woOpenMp);
    free(wOpenMp);
    return 0;
}


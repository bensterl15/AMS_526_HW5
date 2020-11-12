/* Example C code for solving a linear system Ax=b using LAPACK */

#include <stdio.h>

/* Declare function prototype */
extern int sgesv_(int *n, int *nrhs, float *a, int *lda,
                  int *ipiv, float *b, int *ldb, int *info);
/*  -- LAPACK driver routine (version 3.1) --   
    Purpose   
    =======   

    SGESV computes the solution to a real system of linear equations   
       A * X = B,   
    where A is an N-by-N matrix and X and B are N-by-NRHS matrices.   

    The LU decomposition with partial pivoting and row interchanges is   
    used to factor A as   
       A = P * L * U,   
    where P is a permutation matrix, L is unit lower triangular, and U is   
    upper triangular.  The factored form of A is then used to solve the   
    system of equations A * X = B.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the N-by-N coefficient matrix A.   
            On exit, the factors L and U from the factorization   
            A = P*L*U; the unit diagonal elements of L are not stored.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    IPIV    (output) INTEGER array, dimension (N)   
            The pivot indices that define the permutation matrix P;   
            row i of the matrix was interchanged with row IPIV(i).   

    B       (input/output) REAL array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS matrix of right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization   
                  has been completed, but the factor U is exactly   
                  singular, so the solution could not be computed.   

    ===================================================================== */

#define SIZE 3				/* dimension of matrix */

static int solve(float A[][SIZE], float *b) {

    int i,j,info,n,nrhs,lda,ipiv[SIZE],ldb;
    float AT[SIZE][SIZE];
    
    /* Permute matrix */
    for (i=0; i<SIZE; i++) {
        for(j=0; j<SIZE; j++) 
            AT[j][i]=A[i][j];
    }
    
    /* Invoke sgesv_ */
    n = lda = ldb = SIZE; nrhs = 1;
    sgesv_(&n, &nrhs, &AT[0][0], &lda, ipiv, b, &ldb, &info);
    return info;
}

int main(int argc, char **argv)
{
    int i, j, pivot[SIZE], ok;
    float A[SIZE][SIZE], b[SIZE];
    
    /* Matrix A */
    A[0][0]= 1.1;  A[0][1]= 2.2;  A[0][2]=-3.3;
    A[1][0]= 4.4;  A[1][1]=-5.5;  A[1][2]= 6.6;
    A[2][0]=-7.7;  A[2][1]= 8.8;  A[2][2]= 9.9;
    
    /* Define right hand side vector */
    b[0] = 0;
    b[1] = 5.5;
    b[2] = 11;
    
    /* Call warpper function */
    ok = solve(A, b);
    
    /* Print out solution vector x */
    for (j=0; j<SIZE; j++) 
        printf("%g\n", b[j]);
}


/* This is the template file for HW5_6. You need to complete the code by 
 * implementing the parts marked by FIXME. 
 *
 * Submit your completed code to the TA by email. Also submit the plots 
 * and your conclusions of the comparative study.
 *
 * To compile, use the UNIX command 
 *      make
 * in the source directory, and it would then invoke the compilation command 
 * with the included makefile. You may sometimes need to use the command 
 *      make clean
 * before recompiling.
 *
 * This will generate an executable hw5_6. The -g option is optional and 
 * is needed only if you need to debug the program in a debugger such as ddd.
 * The -Wall option would enable compiler's warning messages.
 *
 * Run the program with command
 *      ./hw5_6
 * which would generate a M-file, which you can use to generate the plots by
 * issueing the UNIX command 
 *      make plot
 * in the source directory. It will then generate two PDF files, which you
 * should submit to the TA along with your source code and your conclusions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

/* All type definitions and function prototypes are in this header file. */
#include "hw5_6.h"

#define MAXDEGREE 15
#define MAXPOINTS 30
#define MAXMATRIXSIZE  450
#define MAXRSIZE  225

void print_matrix(Matrix M){
    int m = M.m, n = M.n;
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            printf("%f ", M.vals[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char **argv) 
{   
    /* Buffer vectors */
    double ts[MAXPOINTS];
    double dtemp[MAXDEGREE];

    /* Arrays for storing error and timing results */
    double err_svd[MAXDEGREE-2], err_eig[MAXDEGREE-2];
    double times_svd[MAXDEGREE-2], times_eig[MAXDEGREE-2];

    FILE *fid;
    int m, n, i, j;
    int niter = 1000;

    /* Set n to degree of polynomial */
    for (n=3; n<=MAXDEGREE; ++n) {
        Matrix A; 
        Vector b,x;

        /* Set m to the number of points */
        m = n + n;

        A = allocate_matrix( m, n);
        b = allocate_vector( m);
        x = allocate_vector( n);

        /* Compute Vandermonde matrix A and right-hand side b */
        for (i=0; i<m; i++) {
            ts[i] = (((double)i)) / (((double)(m - 1)));

            A.vals[i][0] = 1.0;
            b.vals[i] = 1.0;
        }

        for (j=1; j<n; j++) {
            for (i=0; i<m; i++) {
                A.vals[i][j] = A.vals[i][j-1] * ts[i];
                b.vals[i] = 1.0 +  b.vals[i] * ts[i];
            }
        }

        tictoc(1);
        /* Invoke least-squares solver using SVD (dgesvd). */
        /* Run for many iterations to obtain more accurate timing results */
        for (j=0; j<niter; ++j) {
            lsqr_dgesvd(x, A, b);
        }
        times_svd[n - 3] =  tictoc(2) / niter;

        for (j=0; j<n; j++) {
            dtemp[j] =  x.vals[j] - 1.0;
        }
        err_svd[n - 3] = norm2_col(dtemp, n);

        tictoc(1);
        /* Invoke least-squares solver using symmetric eigenvalue decomposition 
           of A'*A (dsyev). */
        /* Run for many iterations to obtain more accurate timing results */
        for (j=0; j<niter; ++j) {
            lsqr_dsyev(x, A, b);
        }
        times_eig[n - 3] =  tictoc(2) / niter;

        for (j=0; j<n; j++) {
            dtemp[j] =  x.vals[j] - 1.0;
        }
        err_eig[n - 3] = norm2_col(dtemp, n);

        deallocate_matrix( A);
        deallocate_vector( b);
        deallocate_vector( x);
    }

    /* Write out results into a Matlab file */
    fid = fopen("results.m", "w");
    write_vec(fid, "err_svd", err_svd, MAXDEGREE-2);
    write_vec(fid, "err_eig  ", err_eig, MAXDEGREE-2);
    write_vec(fid, "times_svd", times_svd, MAXDEGREE-2);
    write_vec(fid, "times_eig", times_eig, MAXDEGREE-2);
    fclose(fid);

    /* plotresults; */
    return 0;
}

/*************************************************************
 *
 * FUNCTION: lsqr_dgesvd
 *
 * Least squares via SVD of A.
 *************************************************************/

void  lsqr_dgesvd(
   Vector x, 
   const Matrix A, 
   const Vector b)
{
    int m = A.m;
    int n = A.n;

    Matrix Acopy = allocate_matrix(m,n);
    Matrix U = allocate_matrix(m,n);
    Matrix VT = allocate_matrix(n,n);

    double work[5*MAXPOINTS];
    double S[MAXPOINTS];

    int info;
    int lwork = 5*n;
    int i, j, k;

    /* Make a copy of A to avoid overwriting by dgesvd */
    for (i=0; i<m; ++i) {
        for (j=0; j<n; ++j) {
            Acopy.vals[i][j] = A.vals[i][j];
        }
    }
    
    /* Perform SVD of A.'
       Since A is stored as row-major, it is A' in Fortran convenction. 
       Therefore, we call dgesvd as A'=V S U' to compute SVD, and  
       V is then V' in C convention and U' is U in C convention. */
    dgesvd_("A", "S", &n, &m, &Acopy.vals[0][0], &n, S, &VT.vals[0][0],
            &n, &U.vals[0][0], &n, work, &lwork, &info);

    if (info != 0) {
        printf("ERROR: dgesvd has failed.");
        return;
    }

    /* FIXME: Implement this function */
    /* Compute x = V*inv(S)*U'*b, or equivalent
       x = sum_j ( U(:,j)'*b/sigma_j * V(:,j); */

    /*
    printf("\nA\n");
    print_matrix(A);
    printf("--\n");

    print_matrix(U);

    printf("--S:\n");
    for(i = 0; i < n; i++) printf("%g\n", S[i]);
    printf("--\n");

    printf("V is %d x %d\n", VT.n, VT.m);
    printf("U is %d x %d\n", U.m, U.n);
    printf("A is %d x %d\n", A.m, A.n);
    printf("b is %d\n", b.m);
    printf("x is %d\n", x.m);
    printf("m is %d and n is %d\n", m, n);
    printf("--\n");
    */

    //Allocate vector and initialize to zero:    
    Vector UT_b = allocate_vector(n);
    memset(UT_b.vals, 0, n * sizeof(*UT_b.vals));
    //Now calculate it:
    for(j = 0; j < m; j++){
        for(k = 0; k < n; k++){
            UT_b.vals[k] += U.vals[j][k] * b.vals[j];
        }
    }

    /* Print UT_b vector:
    printf("--\n UT_b");
    for(i = 0; i < n; i++) printf("%g\n", UT_b.vals[i]);
    printf("--\n");
    */

    //Solve S * w = b for w:
    Vector w = allocate_vector(n);
    memset(w.vals, 0, n * sizeof(*w.vals));
    //for(i = 0; i < n; i++) printf("%g\n", w.vals[i]);
    for(i = 0; i < n; i++){
        //if(S[i] < 0.0001) continue;
        //else{
            w.vals[i] = UT_b.vals[i] / S[i];
            //printf("S[i] = %f\n", S[i]);
        //} 
    }
    
    /* Print w:
    printf("--w:\n");
    for(i = 0; i < n; i++) printf("%g\n", w.vals[i]);
    printf("--\n");
    */
    
    // Use memset to initialize x to zero:
    memset(x.vals, 0, x.m * sizeof(*x.vals));
    //Calculate x = Vw here:
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            x.vals[i] += VT.vals[j][i] * w.vals[j];
        }
    }

    /* Print x:
    printf("--x:\n");
    for(i = 0; i < n; i++) printf("%g\n", x.vals[i]);
    printf("--\n");
    */

    deallocate_matrix(U);
    deallocate_matrix(VT);
    deallocate_matrix(Acopy);
    deallocate_vector(UT_b);
    deallocate_vector(w);

    return;
}


/*************************************************************
 *
 * FUNCTION: lsqr_dsyev
 *
 * Least squares via symmetric eigenvalue decomposition of A'*A.
 *************************************************************/

void  lsqr_dsyev(
   Vector x, 
   const Matrix A, 
   const Vector b)
{
    /* FIXME: Implement this function */

    int m = A.m;
    int n = A.n;

    int i, j, k;

    Matrix AT_A = allocate_matrix(n, n);

    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            AT_A.vals[i][j] = 0;
            for(k = 0; k < m; k++){
                AT_A.vals[i][j] += A.vals[k][i] * A.vals[k][j];
            }
        }
    }

    /*
    printf("--AT_A\n");
    print_matrix(AT_A);
    printf("--\n");
    */
    
    Vector AT_b = allocate_vector(n);
    memset(AT_b.vals, 0, n * sizeof(*AT_b.vals));

    for(i = 0; i < n; i++){
        for(j = 0; j < m; j++){
            AT_b.vals[i] += A.vals[j][i] * b.vals[j];
        }
    }

    int lda = n, info, lwork;
    double wkopt;

    //Store the eigenvalues here:
    double work[5*MAXPOINTS];
    double lambda[MAXPOINTS];

    // Solve eigenproblem... Note once again because of Fortran Convention, we get AT_A becomes Q':
    dsyev_( "Vectors", "Upper", &n, &AT_A.vals[0][0], &lda, lambda, work, &lwork, &info );
    // Check for convergence
    if( info > 0 ) {
            printf( "The algorithm failed to compute eigenvalues.\n" );
            exit( 1 );
    }

    //Print lambda:
    /*
    printf("--lambda:\n");
    for(i = 0; i < n; i++) printf("%g\n", lambda[i]);
    printf("--\n");
    
    printf("--\n");
    print_matrix(AT_A);
    printf("--\n");
    */

    //Solve QLQ'x = Ab:
    Vector QTATb = allocate_vector(n);
    memset(QTATb.vals, 0, n * sizeof(*QTATb.vals));

    //LQ'x = Q'(ATb)
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            QTATb.vals[i] += AT_A.vals[i][j] * AT_b.vals[j];
        }
    }

    //Q'x = inv(L)(Q'ATb)
    Vector iLQTATb = allocate_vector(n);
    memset(iLQTATb.vals, 0, n * sizeof(*iLQTATb.vals));
    //Instead of literal diagonal multiplication, just scale each element in the vector:
    for(int i = 0; i < n; i++) iLQTATb.vals[i] = QTATb.vals[i] / lambda[i];

    /* Print x:
    printf("--iLQTATb:\n");
    for(i = 0; i < n; i++) printf("%g\n", iLQTATb.vals[i]);
    printf("--\n");
    */

    memset(x.vals, 0, n * sizeof(*x.vals));
    //x = Q(inv(L)Q'ATb)
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            x.vals[i] += AT_A.vals[j][i] * iLQTATb.vals[j];
        }
    }

    /* Print x:
    printf("--x:\n");
    for(i = 0; i < n; i++) printf("%g\n", x.vals[i]);
    printf("--\n");
    */

    deallocate_matrix(AT_A);
    deallocate_vector(AT_b);
    deallocate_vector(QTATb);

}

/*************************************************************
 *
 * FUNCTION: norm2_col
 *
 * Compute 2-norm of a given vector.
 * Note: It is not a perfect implementation as it does not guard against 
 *       overflow or underflow.
 *************************************************************/

double norm2_col(
    double *vec, 
    int vec_dim)
{
    double v_out;
    int i;

    v_out = 0.0;
    for (i=0; i<vec_dim; i++) {
        v_out = v_out +  vec[i] *  vec[i];
    }
    return sqrt(v_out);
}

/*-----------------------------------------------------------------------------
 * FUNCTION: allocate_matrix - Allocate memory for a given matrix.
 * --------------------------------------------------------------------------*/
/* Disclaimer: The approach used here is not the most efficient way for
 * implementing matrices, but it is adopted here for convenience.
 */
Matrix allocate_matrix(int m, int n) {
    Matrix A;
    double *ptmp;
    int i;

    A.m = m; A.n = n;
    ptmp = (double *)malloc( sizeof(double) * m * n);
    A.vals = (double **)malloc( sizeof(double **) * m);

    for (i=0; i<m; ++i) {
        A.vals[i] = ptmp + i*n;
    }

    return A;
}

/*-----------------------------------------------------------------------------
 * FUNCTION: deallocate_matrix - De-allocate memory for a matrix.
 * --------------------------------------------------------------------------*/
void deallocate_matrix( Matrix A) {
    free( A.vals[0]);
    free( A.vals); A.vals=NULL;
}


/*-----------------------------------------------------------------------------
 * FUNCTION: allocate_vector - Allocate memory for an m vector.
 * --------------------------------------------------------------------------*/
Vector allocate_vector( int m) {
    Vector x;

    x.m = m;
    x.vals = (double *)malloc( sizeof(double) * m);
    return x;
}

/*-----------------------------------------------------------------------------
 * FUNCTION: deallocate_vector - De-allocate the memory for an m vector. 
 * --------------------------------------------------------------------------*/
void deallocate_vector( Vector x) {
    free( x.vals); x.vals=NULL;
}

#include <time.h>
#include <sys/time.h>

/*-----------------------------------------------------------------------------
 * FUNCTION: tictoc - Initialize (1) or get (2) the elapsed time since in seconds
 * --------------------------------------------------------------------------*/
double tictoc (int n)
{
    double    y = 0.0;
    static    struct timeval start_time, end_time; 

    if (n == 1) { 
        gettimeofday(&start_time, NULL);
    }
    else if (n == 2) { 
        gettimeofday(&end_time, NULL);

        y = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec-start_time.tv_usec)*1.e-6;
    }
    return (y); 
}

/*************************************************************
 *
 * FUNCTION: write_vec
 *
 * Write out a vector into file
 *************************************************************/

void  write_vec(FILE *fid, char *name, const double *vec, int n) {
    int i;

    fprintf(fid, "%s=[", name);
    for (i=0; i<n; ++i) {
        fprintf(fid, "%g; ", vec[i]);
    }
    fprintf(fid, "];\n");
}

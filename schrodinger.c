#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>

int main(void) {
    int L = 3, i, j, t;
    int N = 10;
    double M[L][L][N];
    gsl_matrix *E = gsl_matrix_alloc(L, L);
    gsl_vector_complex *eigen = gsl_vector_complex_alloc(L);
    gsl_eigen_nonsymm_workspace *  w = gsl_eigen_nonsymm_alloc(L);

    for(t = 1; t <= N; t++) {
        M[0][0][t-1] = 0;
        M[0][1][t-1] = t;
        M[0][2][t-1] = t;
        M[1][0][t-1] = t;
        M[1][1][t-1] = 0;
        M[1][2][t-1] = 2.0 * t;
        M[2][1][t-1] = t;
        M[2][0][t-1] = t;
        M[2][2][t-1] = 0;

        for(i = 0; i < L; i++) {
            for(j = 0; j < L; j++) {
                gsl_matrix_set(E, i, j, M[i][j][t - 1]);
            }
        }

        gsl_eigen_nonsymm(E, eigen, w); /*diagonalize E which is M at t fixed*/
        printf("#%d\n\n", t);

        for(i = 0; i < L; i++) {
            printf("%d\t%lf\n", i, GSL_REAL(gsl_vector_complex_get(eigen, i)));
        }
        printf("\n");
   }
}
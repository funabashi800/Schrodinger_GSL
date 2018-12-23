#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

const double Q=1.60217733e-19;
const double F=5.0e-2;
const double M=9.1093897e-31;
const double T=300.0;
const double D=1e+16;
const double dx=1e-10;
const int DIM=100;
const double HBAR=1.05457266e-34;
const double KB=1.380658e-23;
const double PI=3.141592654;

double potential(double x){
    return Q*F*x;
}

void Normalise(gsl_vector *Single_eigvector_ptr, FILE *fp_2)
{
    int n;
    double temp;
    double norm_factor, vector_sum=0;
    for (n=0; n<DIM; n++){
        temp = gsl_vector_get(Single_eigvector_ptr, n);
        vector_sum = vector_sum + temp*temp*dx;
    }
    norm_factor = sqrt(vector_sum);
    for (n=0; n<DIM; n++){
        temp = (gsl_vector_get(Single_eigvector_ptr, n))/norm_factor;
        gsl_vector_set(Single_eigvector_ptr, n, temp);
        fprintf(fp_2, "%g, %E\n", dx*n, gsl_vector_get(Single_eigvector_ptr, n)*dx);
        /*output mth value of nth eigenvector*/
    }
}

double subband_contrib(double E_fermi, double E)
{
    double n_sp = 0;
    double DOS = (M*0.17)/(HBAR*HBAR*PI);
    n_sp = DOS*KB*T*log(exp((E_fermi + E)/(KB*T)) + 1);
    return n_sp;
}

void Charge_Density(gsl_matrix *Norm_eigenvector_ptr, gsl_vector *s_eigvector_ptr, gsl_vector *Charge_Density_ptr, double *subband_ptr, double E_fermi, double c_doping_yz, int N, FILE *fp_3)
{
    int n, m;
    double phi, n_sp, temp;
    double c_doping_xyz;
    temp = sqrt(c_doping_yz);
    c_doping_xyz = temp*temp*temp;
    for (n=0; n<N; n++)
    {
        gsl_matrix_get_col(s_eigvector_ptr, Norm_eigenvector_ptr, n);
        n_sp = subband_contrib(E_fermi, subband_ptr[n]);
        printf("ncharge contribution from subband %d is %E coulombs per meter squared\n", n, n_sp);
        for (m=0; m<DIM; m++)
        {
            phi = gsl_vector_get(s_eigvector_ptr, m);
            temp = gsl_vector_get(Charge_Density_ptr, m) - Q*phi*phi*n_sp;
            gsl_vector_set(Charge_Density_ptr, m, temp);
        }
    }
    /*Output*/
    for(m=0; m<DIM; m++){
        //fprintf(fp_3, "n%E", gsl_vector_get(Charge_Density_ptr, m));
        temp = gsl_vector_get(Charge_Density_ptr, m) + Q*c_doping_xyz;
        //temp = temp*dx*dx; /*for Poisson Solver*/
        gsl_vector_set(Charge_Density_ptr, m, temp);
        fprintf(fp_3, "n%E", gsl_vector_get(Charge_Density_ptr, m));
    }
}

int Schrodinger(double E_well, double *subband_ptr,
gsl_eigen_symmv_workspace *Workspace, gsl_matrix *Hamil_ptr, gsl_vector *Eigvalue_ptr, gsl_matrix *Eigvector_ptr,
gsl_vector *s_eigvector_ptr, gsl_vector *Potential_ptr, gsl_matrix *Norm_Eigvector_ptr, FILE *fp_1, FILE *fp_2){
    int n, m, N=0;
    double Element;
    double t_0 = (HBAR*HBAR)/(2.0*dx*dx);
    for(n=0; n<DIM; n++){
        for(m=0; m<DIM; m++){
            if(n==m){
                Element = 2*t_0/(M*0.17)+gsl_vector_get(Potential_ptr, n);
            }
            else if (n == m-1) /*subdiagonal*/
            {
                Element = -1*t_0/(M*0.17);
            }
            else if (n == m+1) /*superdiagonal*/
            {
                Element = -1*t_0/(M*0.17);
            }
            else
            {
                Element = 0;
                /*all other matrix elements*/
            }
            gsl_matrix_set(Hamil_ptr, n, m, Element);
/*transfer value to matrix*/
        }
    }
    gsl_eigen_symmv(Hamil_ptr, Eigvalue_ptr, Eigvector_ptr, Workspace);
    /*solve matrix*/
    gsl_eigen_symmv_sort(Eigvalue_ptr, Eigvector_ptr, GSL_EIGEN_SORT_VAL_ASC);
    /*re-order the results in ascending order*/
    /*output eigenvalues*/
    int count = 0;
    for (n = 0; n < DIM; n++){
        double eigenvalue = gsl_vector_get(Eigvalue_ptr, n);
        if (eigenvalue>E_well) {
            count++;
            fprintf(fp_1, "n%E", (eigenvalue/1.6E-22));
            /*output in meV*/
            subband_ptr[n] = eigenvalue; /*For Fermi-Dirac */
            N = N+1;
            if(count > 3) break;
            /* sum the number of eigenstates*/
        }
    }
    for (n = 0; n < N; n++){  /*retrieve eigenvectors for confined states only*/
        gsl_matrix_get_col(s_eigvector_ptr, Eigvector_ptr, n);
        /*retrieve nth eigenvector*/
        Normalise(s_eigvector_ptr, fp_2);
        /*normalise eigenvector*/
        gsl_matrix_set_col(Norm_Eigvector_ptr, n, s_eigvector_ptr);
    }
    return N;

}

double FD_Bisect(double *Subband_f, double E_fermi_super, double E_fermi_diff, double c_doping_yz, int N)
{
    double n_sp = 0, n_s = 0;
    int n, m;
    /*temperature set to 300K*/
    double DOS = (M*0.17)/(HBAR*HBAR*PI);
    double E_fermi_sub = E_fermi_super + E_fermi_diff;
    double E_fermi_mid = (E_fermi_sub + E_fermi_super)/2;
    for (n=0; n<50; n++){
        for (m = 0; m < N; m++){ /*add contributions from all subbands*/
            n_sp = DOS*KB*T*log(exp((E_fermi_mid - Subband_f[m])/(KB*T)) + 1);
            n_s = n_sp + n_s;
        }
        if (n_s > c_doping_yz){
            E_fermi_super = E_fermi_mid;
            E_fermi_mid = (E_fermi_super + E_fermi_sub)/2;
        }
        else {
            E_fermi_sub = E_fermi_mid;
            E_fermi_mid = (E_fermi_sub + E_fermi_super)/2;
        }
        n_s = 0;
    }
    printf("nPrecise Fermi energy is %E eV\n", E_fermi_mid);
    return E_fermi_mid;
}

double FermiDiracPopulate(double *Subband_f, double c_doping_yz, double E_fermi, int N)
{
    int m;
    double FL_precise, n_sp = 0, n_s = 0;
    /* initial fermi-energy roughly at intrinsic level (@300K)*/
    double d_E_fermi = E_fermi/100;
    double DOS = (M*0.17)/(HBAR*HBAR*PI);
    for (n_s = 0; E_fermi < 0; E_fermi = E_fermi - d_E_fermi){
        printf("n%E\n", E_fermi);
        for (m = 0; m < N; m++){ /*add contributions from all subbands*/
            n_sp = DOS*KB*T*log(exp((E_fermi - Subband_f[m])/(KB*T)) + 1);
            /*NB. subband energies were calculated to be positive (should be negative)*/
            n_s = n_sp + n_s;
        }
        //printf("ncharge density for E_fermi of %E is %E per meter squared", E_fermi, n_s);
        if (n_s > c_doping_yz){
            printf("nRaw Fermi energy is %E eV\n", E_fermi);
            FL_precise = FD_Bisect(Subband_f, E_fermi, d_E_fermi, c_doping_yz, N);
            break;
        }
        n_s = 0;
    }
    return FL_precise;
}

double charge_density_check(gsl_vector *charge_density_ptr)
{
    double charge_tot=0;
    int n;
    for(n=0; n<DIM; n++)
    {
        charge_tot = charge_tot + gsl_vector_get(charge_density_ptr, n);
    }
    return charge_tot;
}

int main(void){
    int n=0, N=0;
    double E_fermi;
    // E_well shows the minimum energy to search eigenvalue
    double c_doping_yz = 1E16, E_well = -1.6022E-20;

    FILE *fp_1, *fp_2;
    FILE *fp_3;
    fp_1=fopen("eigenvalues.txt", "w+");
    fp_2=fopen("eigenvectors.txt", "w+");
    fp_3=fopen("charge_distribution.txt", "w+");
    double *Subband;
    Subband = calloc(DIM, sizeof(double));
    gsl_eigen_symmv_workspace *Workspace = gsl_eigen_symmv_alloc(DIM);
    gsl_matrix *Hamiltonian = gsl_matrix_alloc(DIM, DIM);
    gsl_vector *EigenValue = gsl_vector_alloc(DIM);
    gsl_matrix *EigenVector = gsl_matrix_alloc(DIM, DIM);
    gsl_vector *Potential = gsl_vector_alloc(DIM);
    gsl_matrix *Norm_Eigvector_ptr = gsl_matrix_alloc(DIM, DIM);
    gsl_vector *S_Eigvector_ptr = gsl_vector_alloc(DIM);
    gsl_vector *Charge_Density_ptr = gsl_vector_alloc(DIM);

    // initialize potential
    for (n=0; n<DIM; n++){
        gsl_vector_set(Potential, n, E_well+potential(n));
    }
    N=Schrodinger(E_well, Subband, Workspace, Hamiltonian, EigenValue, EigenVector, S_Eigvector_ptr, Potential, Norm_Eigvector_ptr, fp_1, fp_2);
    E_fermi = FermiDiracPopulate(Subband, c_doping_yz,E_well, N);
    Charge_Density(Norm_Eigvector_ptr, S_Eigvector_ptr, Charge_Density_ptr, Subband, E_fermi, c_doping_yz, N, fp_3);
    double Charge_tot = charge_density_check(Charge_Density_ptr);
    printf("n%E\n", Charge_tot);
    gsl_matrix_free(EigenVector);
    gsl_vector_free(EigenValue);
    gsl_matrix_free(Hamiltonian);
    gsl_eigen_symmv_free(Workspace);
    gsl_matrix_free(Norm_Eigvector_ptr);
    gsl_vector_free(S_Eigvector_ptr);
    free(Subband);
    gsl_vector_free(Charge_Density_ptr);
    return 0;
}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <omp.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_gamma.h>

using namespace std;

const int NF = 11;
const int  Nbins = 100;
const int  Nmass = 100;
const int  Nsigma = 100;

const double  Mmin  =  1.e-2; /*GeV*/
const double  Mmax  =  1.e3;
const double  SigmaMin =   1.e-50; /*cm^-2*/
const double  SigmaMax  =  1.e-40;

const double  Exp = 1e0;

double Nnu[NF][Nbins];
double NChi[Nmass][Nbins];

const double FluxesUn[NF][2] = {{5.25e6, 0.21e6},
                        {7.98e3, 2.39e3},
                        {5.98e10, 0.04e10},
                        {10.5, 2.1},
                        {5.29e6, 1.06e6},
                        {2.78e8, 0.42e8},
                        {2.05e8, 0.35e8},
                        {86., 43.},
                        {4.84e8, 0.15e8},
                        {4.35e9, 0.13e9},
                        {1.44e8, 0.01e8}};

/* B8, Hep, Pp, Atm, F17, N13, 015, DSN, Be7@384.3, Be7@861.3, pep */

double M(int i){
return Mmin*pow(Mmax/Mmin,(double) i/(Nmass-1));
}

double Sigma(int j){
return SigmaMin*pow(SigmaMax/SigmaMin,(double) j/(Nsigma-1));
}

void ReadEvents(const char *filename, double (*array)[Nbins]) {
    FILE *file;
    int count = 0;

    // Open the .dat file for reading
    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        //return EXIT_FAILURE;
    }


    // Read the data
    for (int i = 0; i < NF; i++){
        for (int j = 0; j < Nbins; j++) {
            if (fscanf(file, "%lf", &array[i][j]) != 1) { // Change %lf to %d for integers
                fprintf(stderr, "Error reading data at row %d, column %d\n", i, j);
                fclose(file);
                //return EXIT_FAILURE;
            }
        }
    }
    // Close the file
    fclose(file);

}

void ReadEvents2D(const char *filename, double (*array)[Nbins]){
    FILE *file;

    // Open the file for reading
    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
       //return EXIT_FAILURE;
    }

    // Read the data
    for (int i = 0; i < Nmass; i++) {
        for (int j = 0; j < Nbins; j++) {
            if (fscanf(file, "%lf", &array[i][j]) != 1) { // Change %lf to %d for integers
                fprintf(stderr, "Error reading data at row %d, column %d\n", i, j);
                fclose(file);
                //return EXIT_FAILURE;
            }
        }
    }
    // Close the file
    fclose(file);

}

// Define the function to minimize
double LogLikelihood0(const gsl_vector *v, void *par){
    double x0 = gsl_vector_get(v, 0); //8B flux normalisation
    double x1 = gsl_vector_get(v, 1); //Hep flux normalisation
    double x2 = gsl_vector_get(v, 2); //Atm flux normalisation
    double x3 = gsl_vector_get(v, 3); // flux normalisation
    double x4 = gsl_vector_get(v, 4); // flux normalisation
    double x5 = gsl_vector_get(v, 5); // flux normalisation
    double x6 = gsl_vector_get(v, 6); // flux normalisation
    double x7 = gsl_vector_get(v, 7); // flux normalisation
    double x8 = gsl_vector_get(v, 8); // flux normalisation
    double x9 = gsl_vector_get(v, 9); // flux normalisation
    double x10 = gsl_vector_get(v, 10); // flux normalisation

    double x[NF] = {x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10};

    int *p = (int *) par;

    double B, Bp, S=0, nobs=0, LP0=0;

    for(int i=0; i<Nbins; i++){
        B=0; Bp=0;
        for(int j=0; j<NF; j++){
            B += Exp*x[j]*Nnu[j][i]; // nexp 
            Bp += Exp*FluxesUn[j][0]*Nnu[j][i]; // nobs = nnu + nchi
        }
        
        S = Exp*Sigma(p[1])*NChi[p[0]][i]; //Mass index

        nobs = Bp+S; // ith component

        /*if(B > 0)*/ LP0 += nobs*log(B) - B - gsl_sf_lngamma(nobs+1);
        //else LP0 += 0;
    }

    double SLG = 0;

    for(int i=0; i<NF; i++){
        SLG += -0.5*log(2*M_PI) - log(FluxesUn[i][1]) - 0.5*(x[i]-FluxesUn[i][0])*(x[i]-FluxesUn[i][0])/(FluxesUn[i][1]*FluxesUn[i][1]);
    }

    return -LP0-SLG;
}

double LL1(int i, int j) {
    double B, S, LP1=0;

    for(int k=0; k<Nbins; k++){
        B=0; S=0;
        for(int l=0; l<NF; l++){
            B += Exp*FluxesUn[l][0]*Nnu[l][k];
        }
    
        S = Exp*Sigma(j)*NChi[i][k];

        /*if( S+B > 0)*/ 
        LP1 += (B+S)*log(B+S) - (B+S) - gsl_sf_lngamma(B+S+1);
        /*else LP1 += 0;*/
    }

    double SLG=0; //Calculated using Azimov data set

    for(int l=0; l<NF; l++){
        SLG += -0.5*log(2*M_PI) - log(FluxesUn[l][1]);
    }

    return LP1+SLG;
}

double NWimp(int i, int j) {
    //double B, ntot=0;
    double S=0;

    for(int k=0; k<Nbins; k++){
        //B=0;
        //for(int l=0; l<NF; l++){
            //B += FluxesUn[l][0]*Exp*Nnu[l][k];
        //}
    
        S += Sigma(j)*Exp*NChi[i][k];

    }

    return S;
}

double NTot(int i, int j) {
    double B, S, ntot=0;

    for(int k=0; k<Nbins; k++){
        B=0; S=0;
        for(int l=0; l<NF; l++){
            B += FluxesUn[l][0]*Exp*Nnu[l][k];
        }
        //S = Sigma(j)*Exp*NChi[i][k];

        ntot += S+B;
    }

    return ntot;
}

double NExp(int i, int j, const gsl_vector *v) {
    double x0 = gsl_vector_get(v, 0); //8B flux normalisation
    double x1 = gsl_vector_get(v, 1); //Hep flux normalisation
    double x2 = gsl_vector_get(v, 2); //Atm flux normalisation
    double x3 = gsl_vector_get(v, 3); //DSN flux normalisation
    double x4 = gsl_vector_get(v, 4); // flux normalisation
    double x5 = gsl_vector_get(v, 5); // flux normalisation
    double x6 = gsl_vector_get(v, 6); // flux normalisation
    double x7 = gsl_vector_get(v, 7); // flux normalisation
    double x8 = gsl_vector_get(v, 8); // flux normalisation
    double x9 = gsl_vector_get(v, 9); // flux normalisation
    double x10 = gsl_vector_get(v, 10); // flux normalisation

    double x[NF] = {x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10};

    double B, nexp=0;
    double S=0;

    for(int k=0; k<Nbins; k++){
        B=0;
        for(int l=0; l<NF; l++){
            B += x[l]*Exp*Nnu[l][k];
        }
    
        //S = Sigma(j)*Exp*NChi[i][k];

        nexp += B;
    }

    return nexp;
}


// Clamp values to specified ranges
void Clamp(gsl_vector *v) {
//double Scale = 3;

    if (gsl_vector_get(v, 0) > FluxesUn[0][0]+FluxesUn[0][1]) gsl_vector_set(v, 0, FluxesUn[0][0]+FluxesUn[0][1]);
    if (gsl_vector_get(v, 0) < FluxesUn[0][0]-FluxesUn[0][1]) gsl_vector_set(v, 0, FluxesUn[0][0]-FluxesUn[0][1]);
    if (gsl_vector_get(v, 1) > FluxesUn[1][0]+FluxesUn[1][1]) gsl_vector_set(v, 1, FluxesUn[1][0]+FluxesUn[1][1]);
    if (gsl_vector_get(v, 1) < FluxesUn[1][0]-FluxesUn[1][1]) gsl_vector_set(v, 1, FluxesUn[1][0]-FluxesUn[1][1]);
    if (gsl_vector_get(v, 2) > FluxesUn[2][0]+FluxesUn[2][1]) gsl_vector_set(v, 2, FluxesUn[2][0]+FluxesUn[2][1]);
    if (gsl_vector_get(v, 2) < FluxesUn[2][0]-FluxesUn[2][1]) gsl_vector_set(v, 2, FluxesUn[2][0]-FluxesUn[2][1]);
    if (gsl_vector_get(v, 3) > FluxesUn[3][0]+FluxesUn[3][1]) gsl_vector_set(v, 3, FluxesUn[3][0]+FluxesUn[3][1]);
    if (gsl_vector_get(v, 3) < FluxesUn[3][0]-FluxesUn[3][1]) gsl_vector_set(v, 3, FluxesUn[3][0]-FluxesUn[3][1]);
    if (gsl_vector_get(v, 4) > FluxesUn[4][0]+FluxesUn[4][1]) gsl_vector_set(v, 4, FluxesUn[4][0]+FluxesUn[4][1]);
    if (gsl_vector_get(v, 4) < FluxesUn[4][0]-FluxesUn[4][1]) gsl_vector_set(v, 4, FluxesUn[4][0]-FluxesUn[4][1]);
    if (gsl_vector_get(v, 5) > FluxesUn[5][0]+FluxesUn[5][1]) gsl_vector_set(v, 5, FluxesUn[5][0]+FluxesUn[5][1]);
    if (gsl_vector_get(v, 5) < FluxesUn[5][0]-FluxesUn[5][1]) gsl_vector_set(v, 5, FluxesUn[5][0]-FluxesUn[5][1]);
    if (gsl_vector_get(v, 6) > FluxesUn[6][0]+FluxesUn[6][1]) gsl_vector_set(v, 6, FluxesUn[6][0]+FluxesUn[6][1]);
    if (gsl_vector_get(v, 6) < FluxesUn[6][0]-FluxesUn[6][1]) gsl_vector_set(v, 6, FluxesUn[6][0]-FluxesUn[6][1]);
    if (gsl_vector_get(v, 7) > FluxesUn[7][0]+FluxesUn[7][1]) gsl_vector_set(v, 7, FluxesUn[7][0]+FluxesUn[7][1]);
    if (gsl_vector_get(v, 7) < FluxesUn[7][0]-FluxesUn[7][1]) gsl_vector_set(v, 7, FluxesUn[7][0]-FluxesUn[7][1]);
    if (gsl_vector_get(v, 8) > FluxesUn[8][0]+FluxesUn[8][1]) gsl_vector_set(v, 8, FluxesUn[8][0]+FluxesUn[8][1]);
    if (gsl_vector_get(v, 8) < FluxesUn[8][0]-FluxesUn[8][1]) gsl_vector_set(v, 8, FluxesUn[8][0]-FluxesUn[8][1]);
    if (gsl_vector_get(v, 9) > FluxesUn[9][0]+FluxesUn[9][1]) gsl_vector_set(v, 9, FluxesUn[9][0]+FluxesUn[9][1]);
    if (gsl_vector_get(v, 9) < FluxesUn[9][0]-FluxesUn[9][1]) gsl_vector_set(v, 9, FluxesUn[9][0]-FluxesUn[9][1]);
    if (gsl_vector_get(v, 10) > FluxesUn[10][0]+FluxesUn[10][1]) gsl_vector_set(v, 10, FluxesUn[10][0]+FluxesUn[10][1]);
    if (gsl_vector_get(v, 10) < FluxesUn[10][0]-FluxesUn[10][1]) gsl_vector_set(v, 10, FluxesUn[10][0]-FluxesUn[10][1]);
}

void Maximaser(int i, int j, double *MaxLL0){
    // Set up the minimizer
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;

    int par[2] = {i, j};

    s = gsl_multimin_fminimizer_alloc (T, NF);

    // Initial guess
    x = gsl_vector_alloc (NF);
    gsl_vector_set(x, 0, FluxesUn[0][0]); // Initial x0
    gsl_vector_set(x, 1, FluxesUn[1][0]); // Initial x1
    gsl_vector_set(x, 2, FluxesUn[2][0]); // Initial x2
    gsl_vector_set(x, 3, FluxesUn[3][0]); // Initial x3
    gsl_vector_set(x, 4, FluxesUn[4][0]); // Initial x4
    gsl_vector_set(x, 5, FluxesUn[5][0]); // Initial x5
    gsl_vector_set(x, 6, FluxesUn[6][0]); // Initial x6
    gsl_vector_set(x, 7, FluxesUn[7][0]); // Initial x7
    gsl_vector_set(x, 8, FluxesUn[8][0]); // Initial x8
    gsl_vector_set(x, 9, FluxesUn[9][0]); // Initial x9
    gsl_vector_set(x, 10, FluxesUn[10][0]); // Initial x10

    // Define the function structure
    gsl_multimin_function min_func;
    min_func.n = NF; // Number of variables
    min_func.f = LogLikelihood0;
    min_func.params = par;

    ss = gsl_vector_alloc (NF);
    //gsl_vector_set_all (ss, 1e8);
    gsl_vector_set (ss, 0, 1e6);
    gsl_vector_set (ss, 1, 1e3);
    gsl_vector_set (ss, 2, 1e9);
    gsl_vector_set (ss, 3, 1e0);
    gsl_vector_set (ss, 4, 1e6);
    gsl_vector_set (ss, 5, 1e8);
    gsl_vector_set (ss, 6, 1e8);
    gsl_vector_set (ss, 7, 1e1);
    gsl_vector_set (ss, 8, 1e7);
    gsl_vector_set (ss, 9, 1e8);
    gsl_vector_set (ss, 10, 1e6);


    // Initialize the minimizer
    gsl_multimin_fminimizer_set (s, &min_func, x, ss);

    // Perform the minimization
    size_t iter = 0;
    int status;
    double size;

    do{
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);
        Clamp(s->x); // Clamp values to ranges

        if (status) break; // Error occurred

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (s->size, 1e-4);

        if (status == GSL_SUCCESS){
          printf ("converged to minimum at\n");
        }

    } while (status == GSL_CONTINUE && iter < 1e5);

   // Output the results
    printf("Min at (" );

    for(int k=0; k<NF; k++){
        printf("%g, ", gsl_vector_get(s->x, k));
    }
    printf(")\n");

    printf("Min value Log(L0(x)) = %g\n", -LogLikelihood0(s->x, par));

    //Saving data
    MaxLL0[0] = M(i);
    MaxLL0[1] = Sigma(j);
    MaxLL0[2] = -LogLikelihood0(s->x, par);/*-s->fval;*/
    //for(int k=0; k<NF; k++){
    // MaxLL0[k+5] = gsl_vector_get(s->x, k);} 
    MaxLL0[5] = gsl_vector_get(s->x, 0);
    MaxLL0[6] = gsl_vector_get(s->x, 1);
    MaxLL0[7] = gsl_vector_get(s->x, 2);
    MaxLL0[8] = gsl_vector_get(s->x, 3);
    MaxLL0[9] = gsl_vector_get(s->x, 4);
    MaxLL0[10] = gsl_vector_get(s->x, 5);
    MaxLL0[11] = gsl_vector_get(s->x, 6);
    MaxLL0[12] = gsl_vector_get(s->x, 7);
    MaxLL0[13] = gsl_vector_get(s->x, 8);
    MaxLL0[14] = gsl_vector_get(s->x, 9);
    MaxLL0[15] = gsl_vector_get(s->x, 10);
    //MaxLL0[17] = NExp(i, j, s->x);
    
    // Clean up
    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(x);
    gsl_vector_free(ss);
}

int main() {

//Load Events
    const char *filename = "/Users/mrlamprea/HEPCodes/Coherent/EventsBinned/NuEvents-Ge.dat";
    ReadEvents(filename, Nnu);

    const char *filename2 = "/Users/mrlamprea/HEPCodes/Coherent/EventsBinned/NChiEvents-Ge.dat";
    ReadEvents2D(filename2, NChi);

// Create Grid

double MaxLL0[19];

  FILE *file = fopen("/Users/mrlamprea/HEPCodes/Coherent/DiscoveryLimit/DL-Ge32.dat", "w");
    if (file == NULL) {
        perror("Error opening file");
        return 1; // Exit if file can't be opened
    }

    //#pragma omp parallel for collapse(2) shared(likelihood_grid)
    for (int i = 0; i < Nmass; i++){
        for (int j = 0; j < Nsigma; j++){

            Maximaser(i, j, MaxLL0);

            MaxLL0[3] = LL1(i,j);
            MaxLL0[4] = -2*(MaxLL0[2]-MaxLL0[3]);
            //MaxLL0[16] = NTot(i,j);
            //MaxLL0[18] = NWimp(i,j);

            for (int k = 0; k < NF+4; k++) {
                fprintf(file ,"%e \t", MaxLL0[k]);
            }
        fprintf(file ,"\n");

        }
    }

    fclose(file);
    printf("Data successfully written to file.dat\n");

    return 0;
}
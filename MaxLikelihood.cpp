#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <omp.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_gamma.h>

using namespace std;

#define Nbins   100
#define Nmass   100
#define Nsigma  100

#define Mmin    1.1 /*GeV*/
#define Mmax    1000.
#define SigmaMin    1.e-50 /*cm^-2*/
#define SigmaMax    1.e-40

//double Nnu8B[Nbins];
//double NnuHep[Nbins];
double Nnu[2][Nbins];
double NChi[Nmass][Nbins];
double FluxesUn[2][2] = {{5.58e6, 0.14e6},{8.04e3, 1.30e3}};
/*B8, Hep*/

double M(int i){ 
return Mmin*pow(Mmax/Mmin,(double) i/Nmass);
}

double Sigma(int j){
return SigmaMin*pow(SigmaMax/SigmaMin,(double) j/Nsigma);
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
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < Nbins; j++) {
            if (fscanf(file, "%lf", &array[i][j]) != 1) {
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

    fclose(file);

}

double logGammaN(double *B, double *S){
  if(*B + *S + 1 <= 170){  
        return log(gsl_sf_gamma(*B + *S + 1));
    }
    else  return (*B + *S)*log(*B + *S)-(*B + *S);
}

// Define the function to minimize
double LogLikelihood0(const gsl_vector *v, void *par){
    double x0 = gsl_vector_get(v, 0); //8B flux normalisation
    double x1 = gsl_vector_get(v, 1); //Hep flux normalisation

    double x[2] = {x0, x1};

    int *p = (int *) par;

    double B;    
    double S = 0;
    double LP0 = 0;

    for(int k=0; k<Nbins; k++){
        B=0;
        for(int l=0; l<2; l++){
            B = B + x[l]*Nnu[l][k];
        }
        S = Sigma(p[1])*NChi[p[0]][k]; //Mass index i,
        
        if(S > 0 && B >0) {LP0 = LP0 + (B+S)*log(B)-B-logGammaN(&B,&S);}
        else LP0 = LP0;     
    }
   
    double SLG = 0;

    for(int l=0; l<2; l++){
        SLG = SLG -log(sqrt(2*M_PI)*FluxesUn[l][1]) - 0.5*(x[l]-FluxesUn[l][0])*(x[l]-FluxesUn[l][0])/FluxesUn[l][1]/FluxesUn[l][1];
    }

    //double LG2 = -log(sqrt(2*M_PI)*FluxesUn[1][1]) - 0.5*(y-FluxesUn[1][0])*(y-FluxesUn[1][0])/FluxesUn[1][1]/FluxesUn[1][1];
    //double LL0=-LP0-LG1-LG2;   

    return -LP0-SLG; 
}

double LL1(int i, int j) {
    double B;    
    double S = 0;
    double LP1 = 0;

    for(int k=0; k<Nbins; k++){
        B=0;
        for(int l=0; l<2; l++){
            B = B + FluxesUn[l][0]*Nnu[l][k]; //FluxesUn[l][0]: central values
        }
        S = Sigma(j)*NChi[i][k];
        
        if(S > 0 && B >0){ LP1 = LP1 + (B+S)*log(B+S)-(B+S)-logGammaN(&B,&S);}
        else LP1 = LP1;
    }
   
    //LG[i] = -log(sqrt(2*M_PI)*FluxesUn[i][1]);
    //double LG2 = -log(sqrt(2*M_PI)*FluxesUn[1][1]);

    double SLG=0;
    
    for(int l=0; l<2; l++){
        SLG = SLG -log(sqrt(2*M_PI)*FluxesUn[l][1]);
    }

    return LP1+SLG;
}


// Clamp values to specified ranges
void Clamp(gsl_vector *v) {

    if (gsl_vector_get(v, 0) < 5.44e6) gsl_vector_set(v, 0, 5.44e6);
    if (gsl_vector_get(v, 0) > 5.72e6) gsl_vector_set(v, 0, 5.72e6);
    if (gsl_vector_get(v, 1) < 6.74e3) gsl_vector_set(v, 1, 6.74e3);
    if (gsl_vector_get(v, 1) > 9.34e3) gsl_vector_set(v, 1, 9.34e3);

}

void Maximaser(int i, int j, double *MaxLL0){
    // Set up the minimizer
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;

    int par[2] = {i, j};
 
    s = gsl_multimin_fminimizer_alloc (T, 2);

    // Initial guess
    x = gsl_vector_alloc (2);
    gsl_vector_set(x, 0, FluxesUn[0][0]); // Initial x0: 8B
    gsl_vector_set(x, 1, FluxesUn[1][0]); // Initial x1: Hep
    
    // Define the function structure
    gsl_multimin_function min_func;
    min_func.n = 2; // Number of variables
    min_func.f = LogLikelihood0;
    min_func.params = par;

    ss = gsl_vector_alloc (2);
    gsl_vector_set_all (ss, 10);

    // Initialize the minimizer
    gsl_multimin_fminimizer_set (s, &min_func, x, ss);

    // Perform the minimization
    size_t iter = 0;
    int status;
    //double size;

    do{
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        Clamp(s->x); // Clamp values to ranges

        if (status){ break;} // Error occurred

        //size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size(s->size, 1e-6);

        if (status == GSL_SUCCESS){
            printf ("converged to minimum at\n");
        }

    } while (status == GSL_CONTINUE  && iter < 1000);

    // Output the results
    printf("Minimum at: (x, y) = (%g, %g)\n", gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1));
    printf("Minimum value: LL0(x, y) = %e, %e\n", -s->fval, -LogLikelihood0(s->x, par));

    MaxLL0[0] = M(i);
    MaxLL0[1] = Sigma(j);
    MaxLL0[2] = -s->fval; /*-LogLikelihood0(s->x, par);*/
    MaxLL0[3] = gsl_vector_get(s->x, 0);
    MaxLL0[4] = gsl_vector_get(s->x, 1);
   
    // Clean up
    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(x);
    gsl_vector_free(ss);
}

int main() {

//Load Events 
    const char *filename3 = "/Users/mrlamprea/HEPCodes/Coherent/EventsBinned/NuEvents.dat";
    ReadEvents(filename3, Nnu);

    const char *filename4 = "/Users/mrlamprea/HEPCodes/Coherent/EventsBinned/NChiEvents.dat";
    ReadEvents2D(filename4, NChi);

// Create Grid 

double MaxLL0[7];

  FILE *file = fopen("/Users/mrlamprea/HEPCodes/Coherent/DiscoveryLimit/DiscoveryLimit-50bins-Xe-8B-Hep.dat", "w");
    if (file == NULL) {
        perror("Error opening file");
        return 1; // Exit if file can't be opened
    }

    //#pragma omp parallel for collapse(2) shared(likelihood_grid)
    for (int i = 0; i < Nmass; ++i) {
        for (int j = 0; j < Nsigma; ++j) {

            Maximaser(i, j, MaxLL0); 

            MaxLL0[5] = LL1(i,j);
            MaxLL0[6] = -2*(MaxLL0[2]-MaxLL0[5]);

            for (int i = 0; i < 7; i++) {
                fprintf(file ,"%e \t", MaxLL0[i]);
            }
        fprintf(file ,"\n");

        }
    }

    fclose(file);
    printf("Data successfully written to output.txt\n");

    return 0;
}
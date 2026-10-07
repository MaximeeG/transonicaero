#ifndef WAVESOLVER_H
#define WAVESOLVER_H

//Custom variable that makes setting the correct algorithm easier and code more readable
typedef enum {
    WAVE_BACKWARD,
    WAVE_FORWARD,
    WAVE_LAX,
    WAVE_LAX_WENDROFF,
    WAVE_LEAPFROG,
    WAVE_THETA,
    WAVE_OWN2SPACE4TIME,
    WAVE_OWN4SPACE2TIME
} Algorithm;

// Struct that takes in the general conditions of the simulation
typedef struct {
    unsigned int nx; // number of grid points
    double x0; // boundary value
    double x1; // boundary value
    double c; // constant
    double cfl; // CFL=c*(deltaT/deltaX)
    double theta;
} AlgorithmConfig;

typedef struct {
    AlgorithmConfig *config;

    //step deltas
    double dx;
    double dt;
    double time; //general time variable

    // the following variables are arrays. 
    // C doesnt really differentiate between pointers and arrays which is why they are declared as pointers

    // x vector
    double currentX; // this variable keeps track of the x value of the current state
    double *x;

    double *u_prev; //only for leap frog
    double *u;
    double *u_next;
} WaveSolverState;


// This function initializes the PDE by:
// - copying config data from AlgorithmConfig
// - determining grid spacing
// - alloc memory and build grid
void waveInit(WaveSolverState *wave, AlgorithmConfig *config);

// This function clears all arrays after simulation
void waveClear(WaveSolverState *wave);

// This function inserts the initial condition into the vector u
void waveSetInitialCond(WaveSolverState *wave);

// This function performs one step of the selected algorithm
void waveStep(WaveSolverState *wave, Algorithm algorithm);

// This function writes the current time step into a csv file
void stateWriteToCSV(FILE *outputFile, WaveSolverState *wave);

// This function solves Au=d by using the thomas algorithm.
// - n is the size of the vector u
// - a, b and c are the tridiagonal components of A (n x n)
// - scrath is a temporary storage variable
// d gets overwritten with the solution for u. This saves memory and makes the function a bit more convenient to use in this context
void solveThomas(double a, double b, double c, double *d, double *scratch, unsigned int n);

#endif // WAVESOLVER_H

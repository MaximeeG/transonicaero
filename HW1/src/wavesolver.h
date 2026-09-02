#ifndef WAVESOLVER_H
#define WAVESOLVER_H

//Custom variable that makes setting the correct algorithm easier and code more readable
typedef enum {
    WAVE_BACKWARD,
    WAVE_FORWARD,
    WAVE_LAX,
    WAVE_LAX_WENDROFF,
    WAVE_LEAPFROG,
    WAVE_THETA
} Algorithm;

// Struct that takes in the general conditions of the simulation
typedef struct {
    unsigned int nx; // number of grid points
    double x0; // boundary value
    double x1; // boundary value
    double c; // constant
    double cfl; // CFL=c*(deltaT/deltaX)
    double theta;// see HW1 hybrid explicit-implicit
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
    double *x;

    double *u_prev; //only for leap frog
    double *u;
    double *u_next;
} WaveSolverState;

// All struct arguments are passed by reference to improve efficiency:

// This function initializes the PDE by:
// - copying config data from AlgorithmConfig
// - determining grid spacing
// - alloc memory and build grid
void waveInit(WaveSolverState *wave, AlgorithmConfig *config);

// This function...
void waveClear(WaveSolverState *wave);

// This function inserts the initial condition into the vector u
void waveSetInitialCond(WaveSolverState *wave, AlgorithmConfig *config);

// This function...
void waveStep(WaveSolverState *wave, Algorithm *algorithm);

// This function...
void waveWriteToCSV(void);

#endif // WAVESOLVER_H

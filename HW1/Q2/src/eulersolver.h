#ifndef EULERSOLVER_H
#define EULERSOLVER_H

#include <stdio.h> 

typedef enum {
    MACCORMACK,
    BEAM_WARMING
} EulerAlgorithm;

typedef enum {
    SUPERSONIC_OUTLET,
    SUBSONIC_OUTLET
} OutletType;

typedef struct {
    unsigned int nx; 
    double x0;      // 0
    double x1;      // 10
    double cfl;     // Target CFL = max(|u| + a) * dt / dx

    double gamma;   // Specific heat ratio 
    double rho_in;  
    double p_in;    
    double mach_in; // 1.25

    OutletType outlet_type;
    double back_pressure_ratio; // p_out / p_in: 1.9 for the subsonic case;

    EulerAlgorithm algorithm;
    unsigned int max_steps;     // Maximum number of time-marching iterations
    double residual_tolerance; // Threshold for the normalized steady residual

} AlgorithmConfig;

typedef struct {
    AlgorithmConfig *config;

    double dx;
    double dt;
    double time;

    double *x;
    double *area;       // A(x)
    double *darea_dx;   // A'(x)

    // Conserved state: each array contains config->nx rows of 3 doubles.
    // Q[i][0] = rho * A       (mass per unit length)
    // Q[i][1] = rho * u * A   (momentum per unit length)
    // Q[i][2] = e * A         (total energy per unit length)
    // Here e is total energy per unit volume, and u is fluid velocity.
    double (*Q)[3];
    double (*Q_predictor)[3]; // MacCormack predictor stage
    double (*Q_next)[3];      // Updated state for the next time step

} EulerSolverState;


// Function that initializes arrays and structs with the correct 
// data to start the simulation
void eulerInit(EulerSolverState *state, AlgorithmConfig *config);


// Function that clears all arrays after simulation
void eulerClear(EulerSolverState *state);


// Set state to the initial condition before starting the simulation
void eulerSetInitialCond(EulerSolverState *state);


// Performs one simulation step with either MacCormack or Beam-Warming algorithm
void eulerStep(EulerSolverState *state);


// Writes data to a csv for analysis
void stateWriteToCSV(FILE *outputFile, const EulerSolverState *state);
 

#endif

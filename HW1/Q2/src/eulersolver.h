#ifndef EULERSOLVER_H
#define EULERSOLVER_H

// Struct that takes in the general conditions of the simulation
typedef struct {
    unsigned int nx; // number of grid points
    double x0; // boundary value
    double x1; // boundary value
    double c; // constant
    double cfl; // CFL=c*(deltaT/deltaX)
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
} EulerSolverState;

void eulerInit(void);

void eulerSetInitialCond(void);

void eulerStep(void);
 
void stateWriteToCSV(FILE *outputFile, EulerSolverState *state, AlgorithmConfig *config);


#endif
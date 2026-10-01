#include "eulersolver.h"

void waveInit(EulerSolverState *state, AlgorithmConfig *config){
    // (x1-x0)/(nx-1) to calculate step size
    // save step size to current state struct

    state->dx = ((config->x1) - (config->x0)) / (config->nx - 1);
    state->dt = (config->cfl * state->dx) / config->c;

    state->x      = malloc(config->nx * sizeof(double));
    // calloc() is basically the same as malloc() except it initilizes all values as 0
    state->u_prev = calloc(config->nx, sizeof(double));
    state->u      = calloc(config->nx, sizeof(double));
    state->u_next = calloc(config->nx, sizeof(double));
    
    //loop that fills x vector
    for(int i = 0; i < config->nx; i++){ // the < is really important here to avoid seg faults
    
        state->x[i] = config->x0 + (i * state->dx);
    }

    // set time to zero
    state->time = 0.0;
    
    // after all the vectors are initialized, the initial conditions are entered
    eulerSetInitialCond(state, config);
}

void eulerSetInitialCond(EulerSolverState *state, AlgorithmConfig *config){
    
}
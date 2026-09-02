#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "wavesolver.h"

#define PI 3.14159265358979323846
#define E 2.71828182845904523536

void waveInit(WaveSolverState *wave, AlgorithmConfig *config){
    // (x1-x0)/(nx-1) to calculate step size
    // save step size to current state struct
    wave->dx = ((config->x1) - (config->x0)) / (config->nx - 1);

    // allocate memory. Why? Because memory needs to be allocated outside of this function.
    // if we don't malloc at the location of the referenced struct (config in this case),
    // the data calculated in this function will cease to exist after waveInit() returns.
    wave->x      = malloc(config->nx * sizeof(double));
    // calloc() is basically the same as malloc() except it initilizes all values as 0
    wave->u_prev = calloc(config->nx, sizeof(double));
    wave->u      = calloc(config->nx, sizeof(double));
    wave->u_next = calloc(config->nx, sizeof(double));
    
    //loop that fills x vector
    for(int i = 0; i < config->nx; i++){ // the < is really important here to avoid seg faults
        // x[i] = x0[i] * dx but written in weird C syntax
        wave->x[i] = (config->x0 + i) * wave->dx;
        //printf("%lf ", vecSizeNX[i]); // for debug
    }

    // set time to zero
    wave->time = 0.0;
    
    waveSetInitialCond(wave, config);

}
void waveClear(WaveSolverState *wave){

}

void waveSetInitialCond(WaveSolverState *wave, AlgorithmConfig *config){
    // The initial condition is hard coded in this function
    // Changing the initial condition therefore requires changing this function
    // Initial condition in HW1:
    // u=1 for 0.5 ≤ x ≤ 1
    // u=0 on all other points

    // this loops through the entire array. 
    // its probably not the most efficient way to do this but I'm too lazy to think of a better option now
    for(int i = 0; i < config->nx; i++){
        if (wave->x[i] >= 0.5 && wave->x[i] <= 1.0)
        {
            wave->u[i] = 1.0;
        }
    }

}
void waveStep(WaveSolverState *wave, Algorithm *algorithm){

}

void waveWriteToCSV(void);
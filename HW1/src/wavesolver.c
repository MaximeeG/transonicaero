#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "wavesolver.h"

#define PI 3.14159265358979323846
#define E 2.71828182845904523536

void waveInit(WaveSolverState *wave, AlgorithmConfig *config){
    // (x1-x0)/(nx-1) to calculate step size
    double stepSize = ((config->x1) - (config->x0)) / (config->nx - 1);
    
    
    //loop that fills x vector
    double vecSizeNX[config->nx];
    for(int i = 0; i < config->nx; i++){ // the < is really important here to avoid seg faults
        vecSizeNX[i] = stepSize*i;
        //printf("%lf ", vecSizeNX[i]); // for debug
    }
    wave->x = vecSizeNX; //save xVec to the state struct


    // initialize all u vectors with zeroes
    for(int i = 0; i < config->nx; i++){
        vecSizeNX[i] = 0;
        //printf("%lf ", vecSizeNX[i]); // for debug
    }
    wave->u_prev = vecSizeNX;
    wave->u = vecSizeNX;
    wave->u_next = vecSizeNX;
    
    // The beginning of time...
    wave->time = 0;

    


}
void waveClear(WaveSolverState *wave){

}

void waveSetInitialCond(WaveSolverState *wave){

}
void waveStep(WaveSolverState *wave, Algorithm *algorithm){

}

void waveWriteToCSV(void);
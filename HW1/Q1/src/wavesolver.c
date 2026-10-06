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
    wave->dt = (config->cfl * wave->dx) / config->c;

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
        // JAKOB: Changed next line. Same result for x0 = 0, but otherwise it would be wrong.
        wave->x[i] = config->x0 + (i * wave->dx);
        //printf("%lf ", vecSizeNX[i]); // for debug
    }

    // set time to zero
    wave->time = 0.0;
    
    // after all the vectors are initialized, the initial conditions are entered
    waveSetInitialCond(wave, config);
}

void waveClear(WaveSolverState *wave){
    // not sure how this will be needed later so it is emtpy for now
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
        } else {
            wave->u[i] = 0.0;
        }
    }
}

void solveThomas(double a, double b, double c, double *d, double *scratch, unsigned int n){

    if (n == 0) {
        return;
    }

    // first row
    scratch[0] = c / b;
    d[0] /= b;

    // forward elimination
    for (unsigned int i = 1; i < n; i++) {
        double denominator = b - a * scratch[i - 1];

        scratch[i] = c / denominator;
        d[i] = (d[i] - a * d[i - 1]) / denominator;
    }

    // back substitution
    for (unsigned int i = n - 1; i > 0; i--) {
        d[i - 1] -= scratch[i - 1] * d[i];
    }
}

void waveStep(WaveSolverState *wave, AlgorithmConfig *config, Algorithm algorithm){

    int i;
    double *temp;

    switch (algorithm)
    {
    case WAVE_BACKWARD:

        // perform algorithm
        for(i = 1; i < config->nx; i++){
            wave->u_next[i] = wave->u[i] - config->cfl * (wave->u[i] - wave->u[i-1]);
        }

        // increment time and space variables
        wave->time += wave->dt;
        wave->currentX += wave->dx;

        // swap u vectors
        temp = wave->u;
        wave->u = wave->u_next;
        wave->u_next = temp;


        break;
    case WAVE_FORWARD:

        // perform algorithm
        // Changed last index to nx - 1
        for(i = 1; i < config->nx - 1; i++){
            wave->u_next[i] = wave->u[i] - config->cfl * (wave->u[i+1] - wave->u[i]);
        }

        // increment time and space variables
        wave->time += wave->dt;
        wave->currentX += wave->dx;

        // swap u vectors
        temp = wave->u;
        wave->u = wave->u_next;
        wave->u_next = temp;

        break;

    case WAVE_LAX_WENDROFF:

        for (i = 1; i < config->nx - 1; i++) {
            wave->u_next[i] =
            wave->u[i]
            - 0.5 * config->cfl
            * (wave->u[i+1] - wave->u[i-1])
            + 0.5 * config->cfl * config->cfl
            * (wave->u[i+1] - 2.0 * wave->u[i] + wave->u[i-1]);
        }

        // increment time and space variables
        wave->time += wave->dt;
        wave->currentX += wave->dx;

        // swap u vectors
        temp = wave->u;
        wave->u = wave->u_next;
        wave->u_next = temp;

        break;

    case WAVE_LAX:
        for (i = 1; i < config->nx - 1; i++) {

        wave->u_next[i] =
        0.5 * (wave->u[i+1] + wave->u[i-1])
        - 0.5 * config->cfl
        * (wave->u[i+1] - wave->u[i-1]);
        }
        
        // increment time and space variables
        wave->time += wave->dt;
        wave->currentX += wave->dx;

        // swap u vectors
        temp = wave->u;
        wave->u = wave->u_next;
        wave->u_next = temp;

        break;

    case WAVE_LEAPFROG:

        //use Lax-Wendroff to initialize since there is no previous time level
        if (wave->time == 0.0) {
            for (i = 1; i < config->nx - 1; i++) {
                wave->u_next[i] =
                wave->u[i]
                - 0.5 * config->cfl
                * (wave->u[i+1] - wave->u[i-1])
                + 0.5 * config->cfl * config->cfl
                * (wave->u[i+1] - 2.0 * wave->u[i] + wave->u[i-1]);
            }
        } else {
            // leapfrog
            for (i = 1; i < config->nx - 1; i++) {
                wave->u_next[i] =
                wave->u_prev[i]
                - config->cfl
                * (wave->u[i+1] - wave->u[i-1]);
            }
        }

        // increment time and space variables
        wave->time += wave->dt;
        wave->currentX += wave->dx;

        // swap u vectors
        temp = wave->u_prev;
        wave->u_prev = wave->u;
        wave->u = wave->u_next;
        wave->u_next = temp;

        break;

    case WAVE_THETA:
    { // this case is in its own scope to avoid VLA errors. 

        double theta = config->theta;
        double d[config->nx];
        double scratch[config->nx];

        
        double a = (-1) * theta * (config->cfl/2.0);
        double b = 1.0; 
        double c = theta * (config->cfl/2.0);

        for(int i = 0; i < config->nx - 1; i++){
            //calculate RHS
            d[i] = wave->u[i] - (1.0 - theta) * (config->cfl / 2.0) * (wave->u[i+1] - wave->u[i-1]); 
        }

        // copy the solution and set zero boundary values
        wave->u_next[0] = 0.0;
        wave->u_next[config->nx - 1] = 0.0;


        solveThomas(a, b, c, &d[1], scratch, config->nx - 2);

        for (i = 1; i < config->nx - 1; i++) {
            wave->u_next[i] = d[i];
        }

        // increment time and space variables
        wave->time += wave->dt;
        wave->currentX += wave->dx;

        // swap u vectors
        temp = wave->u;
        wave->u = wave->u_next;
        wave->u_next = temp;

        break;
    }

case WAVE_OWN2SPACE4TIME:
    {
        double cfl = config->cfl;
        double cfl2 = cfl * cfl;
        double cfl3 = cfl2 * cfl;
        double cfl4 = cfl3 * cfl;

        for (i = 2; i < config->nx - 2; i++) {
            wave->u_next[i] = wave->u[i]
                - (cfl / 2.0) * (wave->u[i+1] - wave->u[i-1])
                + (cfl2 / 2.0) * (wave->u[i+1] - 2.0 * wave->u[i] + wave->u[i-1])
                - (cfl3 / 12.0) * (wave->u[i+2] - 2.0 * wave->u[i+1] + 2.0 * wave->u[i-1] - wave->u[i-2])
                + (cfl4 / 24.0) * (wave->u[i+2] - 4.0 * wave->u[i+1] + 6.0 * wave->u[i] - 4.0 * wave->u[i-1] + wave->u[i-2]);
        }
        
        // Handle boundaries
        wave->u_next[0] = 0.0;
        wave->u_next[1] = 0.0;
        wave->u_next[config->nx - 2] = 0.0;
        wave->u_next[config->nx - 1] = 0.0;

        // increment time and space variables
        wave->time += wave->dt;
        wave->currentX += wave->dx;

        // swap u vectors
        temp = wave->u;
        wave->u = wave->u_next;
        wave->u_next = temp;

        break;
    }

    case WAVE_OWN4SPACE2TIME:
    {
        double cfl = config->cfl;
        double cfl2 = cfl * cfl;

        for (i = 2; i < config->nx - 2; i++) {
            wave->u_next[i] = wave->u[i]
                - (cfl / 12.0) * (-wave->u[i+2] + 8.0 * wave->u[i+1] - 8.0 * wave->u[i-1] + wave->u[i-2])
                + (cfl2 / 24.0) * (-wave->u[i+2] + 16.0 * wave->u[i+1] - 30.0 * wave->u[i] + 16.0 * wave->u[i-1] - wave->u[i-2]);
        }

        // Handle boundaries 
        wave->u_next[0] = 0.0;
        wave->u_next[1] = 0.0;
        wave->u_next[config->nx - 2] = 0.0;
        wave->u_next[config->nx - 1] = 0.0;

        // increment time and space variables
        wave->time += wave->dt;
        wave->currentX += wave->dx;

        // swap u vectors
        temp = wave->u;
        wave->u = wave->u_next;
        wave->u_next = temp;

        break;
    }
    
    default:
        break;
    }
}

void stateWriteToCSV(FILE *outputFile, WaveSolverState *wave, AlgorithmConfig *config){
    // structure of the csv file:
    // header: t,c,u (the c value is constant for all time steps. it is only exported to make plotting the analytical solution easier)

    fprintf(outputFile, "%lf,%lf,", wave->time, config->c);
    
    fprintf(outputFile, "[");
    // the loop leaves out the last element
    for (int i = 0; i < (config->nx) - 2; i++)
    {
        fprintf(outputFile, "%lf,", wave->u[i]);
    }
    // last element added here with brackets and line break
    fprintf(outputFile, "%lf]\n", wave->u[config->nx - 1]);
    

}
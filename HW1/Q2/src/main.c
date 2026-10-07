#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eulersolver.h"


int main(void){

    EulerAlgorithm algorithm = BEAM_WARMING;
    OutletType outlet = SUBSONIC_OUTLET;

    // filename setup here:
    char algName[20] = "BEAM_WARMING";
    char outletName[10] = "SUB";
    double cfl = 0.5;
    long unsigned int nx = 401; // must be greater than 5

    // assemble file name and open file
    char outputFileName[100];

    snprintf(
        outputFileName,
        sizeof(outputFileName),
        "plots/data/%s_%s_CFL_%d_NX_%lu.csv",
        algName,
        outletName,
        (int)(cfl*100), // CFL in filename ist factored x100 so that there is no dot
        nx 
    );

    FILE *outputFile = fopen(outputFileName, "w");

    if(outputFile == NULL){
        printf("main.c: output file could not be created.");
        return 1;
    }

    AlgorithmConfig config = {
        .nx = nx,
        .x0 = 0.0,
        .x1 = 10.0,
        .cfl = cfl,

        // nondimensional inlet values, assuming gamma = 1.4
        .gamma = 1.4,
        .rho_in = 1.0,
        .p_in = 1.0,
        .mach_in = 1.25,

        .outlet_type = outlet,
        .back_pressure_ratio = 1.9, // only used for SUBSONIC_OUTLET

        .algorithm = algorithm,
        .max_steps = 30000,
        .residual_tolerance = 1e-8
    };

    // initializes struct with zeroes/NULL pointers
    EulerSolverState state = {0};

    eulerInit(&state, &config);
    stateWriteToCSV(outputFile, &state);

    double residual = INFINITY;
    unsigned int step = 0;

    while (step < config.max_steps && residual > config.residual_tolerance) {

        eulerStep(&state);
        stateWriteToCSV(outputFile, &state);
        residual = 0.0;

        for (unsigned int i = 1; i < config.nx - 1; i++) {
            for (int k = 0; k < 3; k++) {
                double change = fabs(state.Q[i][k] - state.Q_next[i][k]);
                double scale = fabs(state.Q[0][k]);
                double value = change / (state.dt * scale);

                residual = fmax(residual, value);
            }
        }
        step++;
    }

    if (residual <= config.residual_tolerance) {
        printf("Converged after %u steps.\n", step);
    } else {
        printf("Maximum steps reached without convergence.\n");
    }

    fclose(outputFile);
    eulerClear(&state);

    return 0;
}
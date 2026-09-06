#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wavesolver.h"


static Algorithm parseAlgorithm(const char *name)
{
    if (strcmp(name, "BACKWARD") == 0)
        return WAVE_BACKWARD;

    if (strcmp(name, "FORWARD") == 0)
        return WAVE_FORWARD;

    if (strcmp(name, "LAX") == 0)
        return WAVE_LAX;

    if (strcmp(name, "LAX_WENDROFF") == 0)
        return WAVE_LAX_WENDROFF;

    if (strcmp(name, "LEAPFROG") == 0)
        return WAVE_LEAPFROG;

    // invalid value
    return -1;
}

int main(int argc, char *argv[])
{

    // check if correct number of arguments was passed
    if (argc != 4) {
        printf("main.c: Argument mismatch. Usage: %s <algorithm> <cfl> <nx>\n", argv[0]);
        return 1;
    }

    // get selected algorithm from arguments
    Algorithm activeAlg = parseAlgorithm(argv[1]);
    // get selected CFL and nx from arguments
    double cfl = strtod(argv[2], NULL);
    unsigned long nx = strtoul(argv[3], NULL, 10);

    // integrity check
    if (activeAlg == -1) {
        printf("main.c: Unknown algorithm: %s\n", argv[1]);
        return 1;
    }
    if(nx < 2){
        printf("main.c: nx cannot be lower than 2.\n");
        return 1;
    }

    
    // initialize structs with placeholder data
    AlgorithmConfig config = {
        .nx = nx,
        .x0 = 0.0,
        .x1 = 3.14159265,
        .c = 1.0,
        .cfl = cfl,
        .theta = 0.0
    };
    WaveSolverState state = {
        &config,
        0,
        0,
        0,
        0,
        0,
        0,
        0
    };

    // assemble file name and open file
    char outputFileName[100];

    snprintf(
        outputFileName,
        sizeof(outputFileName),
        "HW1/plots/data/WAVE_%s_CFL_%d_NX_%lu.csv",
        argv[1],
        (int)(cfl*100), // CFL in filename ist factored x100 so that there is no dot
        nx 
    );

    FILE *outputFile = fopen(outputFileName, "w");

    // integrity check
    if(outputFile == NULL){
        printf("main.c: Data output file creation failed.");
        return 1;
    }

    fprintf(outputFile, "t,c,u\n");

    waveInit(&state, &config);
    // write t=0 state to file
    stateWriteToCSV(outputFile, &state, &config);
    
    // stop condition: mid wave reaches x=2.5
    // calculate corresponding time:
    double targetTime = (2.5 - 0.75) / config.c;

    do {
        waveStep(&state, &config, activeAlg);
        stateWriteToCSV(outputFile, &state, &config);
    } while (state.time < targetTime);
    
    fclose(outputFile);
    waveClear(&state); //function not implented yet (TO DO)

    return 0;
}
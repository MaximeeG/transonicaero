#include <stdio.h>
#include <string.h>
#include "gnuplot.h"
#include "wavesolver.h"

#define TO_STRING(variable) #variable

int main(void)
{
    // // File decleration
    // FILE *gp = create_gnuplot();
    // FILE *datafile = fopen("build/data.csv", "w");

    // // File integrity check
    // if (gp == NULL || datafile == NULL) {
    //     printf("File creation failed.");
    //     return 1;
    // }

    // // Data file is edited here 
    // fprintf(datafile, "x,y\n"); // csv header
    
    // // This is where the magic happens
    // int y = 0;
    // for(int x = 0; x < 20; x++){
    //     y = 0.5*x*x*x;
    //     fprintf(datafile, "%d, %d\n", x, y);

    // }

    // // Save plot
    // save_plot(gp, "build/data.csv", "build/plot.png");

    // // Close files, avoid memory leak and stuff
    // fclose(datafile);
    // close_gnuplot(gp);

    AlgorithmConfig config = {
        100,
        0,
        3.14159265,
        1, // c
        0.6, //CFL
        0
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

    // I could maybe write a function that manages all this file stuff
    Algorithm activeAlg = WAVE_BACKWARD;
    FILE *outputFile = fopen("HW1/plots/WAVE_BACKWARD.csv", "w");

    if(outputFile == NULL){
        printf("Output file creation failed.");
        return 1;
    }

    fprintf(outputFile, "t,u\n");

    waveInit(&state, &config);
    
    
    // stop condition: mid wave reaches x=2.5
    // calculate corresponding time:
    double targetTime = (2.5 - 0.75) / config.c;

    do {
        waveStep(&state, &config, activeAlg);
        stateWriteToCSV(outputFile, &state, &config, activeAlg);
    } while (state.time < targetTime);
    
    
    

    fclose(outputFile);

    return 0;
}
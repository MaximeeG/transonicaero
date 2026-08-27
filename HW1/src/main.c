#include <stdio.h>
#include "gnuplot.h"

int main(void)
{
    // File decleration
    FILE *gp = create_gnuplot();
    FILE *datafile = fopen("build/data.csv", "w");

    // File integrity check
    if (gp == NULL || datafile == NULL) {
        printf("File creation failed.");
        return 1;
    }

    // Data file is edited here 
    fprintf(datafile, "x,y\n"); // csv header
    
    // This is where the magic happens
    int y = 0;
    for(int x = 0; x < 20; x++){
        y = 0.5*x*x*x;
        fprintf(datafile, "%d, %d\n", x, y);

    }

    // Save plot
    save_plot(gp, "build/data.csv", "build/plot.png");

    // Close files, avoid memory leak and stuff
    fclose(datafile);
    close_gnuplot(gp);

    return 0;
}
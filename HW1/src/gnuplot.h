#ifndef GNUPLOT_H
#define GNUPLOT_H

#include <stdio.h>

/*Creates new gnuplot instance.
* Returns FILE pointer to gnuplot instance or NULL if error occurred.
*/
FILE *create_gnuplot(void);

/*Simply saves the plot to a file
* gets pointer to gnuplot instance, path to data file (plot data) and path to output file (plot image) as arguments
*/
void save_plot(
    FILE *gnuplot,
    const char *data_file,
    const char *output_file
);

/*Closes gnuplot existing passed as argument*/
void close_gnuplot(FILE *gnuplot);

#endif
#include "gnuplot.h"

FILE *create_gnuplot(void)
{
    FILE *gnuplot = popen("gnuplot -persist", "w");

    if (gnuplot == NULL) {
        fprintf(stderr, "Error: could not start gnuplot.\n");
        return NULL;
    }

    return gnuplot;
}

void save_plot(FILE *gnuplot, const char *data_file, const char *output_file)
{
    if (gnuplot == NULL) {
        return;
    }

    fprintf(gnuplot, "set terminal pngcairo size 1000,600\n");
    fprintf(gnuplot, "set output '%s'\n", output_file);

    fprintf(gnuplot, "set title 'My Plot'\n");
    fprintf(gnuplot, "set xlabel 'x'\n");
    fprintf(gnuplot, "set ylabel 'y'\n");

    fprintf(
        gnuplot,
        "plot '%s' using 1:2 every ::1 with linespoints title 'Data'\n", /*plots content from csv file, skips first header line*/
        data_file
    );

    /* Important: finish writing the image file. */
    fprintf(gnuplot, "unset output\n");

    fflush(gnuplot);
}

void close_gnuplot(FILE *gnuplot)
{
    if (gnuplot != NULL) {
        pclose(gnuplot);
    }
}
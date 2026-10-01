#ifndef EULERSOLVER_H
#define EULERSOLVER_H

typedef enum {
    MACCORMACK,
    BEAM_WARMING
} EulerAlgorithm;

typedef enum {
    SUPERSONIC_OUTLET,
    SUBSONIC_OUTLET
} OutletType;

// Struct that takes in the general conditions of the simulation
typedef struct {
    unsigned int nx; // Number of grid points (at least 2)
    double x0;      // Domain start: 0 for this homework
    double x1;      // Domain end: 10 for this homework
    double cfl;     // Target CFL = max(|u| + a) * dt / dx

    double gamma;   // Specific heat ratio (greater than 1)
    double rho_in;  // Positive inlet density
    double p_in;    // Positive inlet static pressure
    double mach_in; // Inlet Mach number: 1.25 for this homework

    OutletType outlet_type;
    double back_pressure_ratio; // p_out / p_in: 1.9 for the subsonic case;
                                // unused for a supersonic outlet

    EulerAlgorithm algorithm;
    unsigned int max_steps;     // Maximum number of time-marching iterations
    double residual_tolerance; // Threshold for the normalized steady residual
} AlgorithmConfig;

typedef struct {
    AlgorithmConfig *config;

    // Grid spacing, time step, and current simulation time.
    double dx;
    double dt;
    double time;

    // Geometry: each array contains config->nx values.
    double *x;
    double *area;       // A(x)
    double *darea_dx;   // dA/dx

    // Conserved state: each array contains config->nx rows of 3 doubles.
    // Q[i][0] = rho * A       (mass per unit length)
    // Q[i][1] = rho * u * A   (momentum per unit length)
    // Q[i][2] = e * A         (total energy per unit length)
    // Here e is total energy per unit volume, and u is fluid velocity.
    double (*Q)[3];
    double (*Q_predictor)[3]; // MacCormack predictor stage
    double (*Q_next)[3];      // Updated state for the next time step

    // Recover primitive variables when needed rather than storing copies:
    // rho = Q[i][0] / area[i]
    // u   = Q[i][1] / Q[i][0]
    // e   = Q[i][2] / area[i]
    // p   = (gamma - 1) * (e - 0.5 * rho * u * u)
} EulerSolverState;

void eulerInit(EulerSolverState *state, AlgorithmConfig *config);

void eulerSetInitialCond(EulerSolverState *state);

void eulerStep(EulerSolverState *state);

void stateWriteToCSV(FILE *outputFile, const EulerSolverState *state);
 

#endif

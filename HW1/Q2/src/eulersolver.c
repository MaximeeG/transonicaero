#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "eulersolver.h"

// artificial dissipation settings. The pressure sensor turns on second-order
// differences around shocks. Fourth-order differences damp grid oscillations
// in the smooth parts without adding first-order spatial error there.
#define DISSIPATION_2 1.0
#define DISSIPATION_4 0.02

// where does this come from? Does it really have to be its own function?
static double pressure(const double Q[3], double area, double gamma){
    // Q contains area-weighted values, so divide by A to get pressure
    return (gamma - 1.0) * (Q[2] - 0.5 * Q[1] * Q[1] / Q[0]) / area;
}


static void setBoundaryCond(const EulerSolverState *state, double (*Q)[3]){
    AlgorithmConfig *config = state->config;
    unsigned int last = config->nx - 1;

    // supersonic inlet: all three quantities are imposed
    double u = config->mach_in * sqrt(config->gamma * config->p_in / config->rho_in);
    double e = config->p_in / (config->gamma - 1.0) + 0.5 * config->rho_in * u * u;

    Q[0][0] = config->rho_in * state->area[0];
    Q[0][1] = config->rho_in * u * state->area[0];
    Q[0][2] = e * state->area[0];

    // extrapolate the primitive state, not Q itself. The area is different
    // at the last two points, even if rho, u and p are the same.
    double areaRatio = state->area[last] / state->area[last - 1];

    for (int k = 0; k < 3; k++) {
        Q[last][k] = areaRatio * Q[last - 1][k];
    }

    if (config->outlet_type == SUBSONIC_OUTLET) {
        // keep extrapolated density and velocity, but impose outlet pressure
        double p_out = config->back_pressure_ratio * config->p_in;

        // this here needs more documentation but I am too lazy rn
        Q[last][2] = p_out * state->area[last] / (config->gamma - 1.0) + 0.5 * Q[last][1] * Q[last][1] / Q[last][0];
    }
}

// what does this function do? 
// I think it only calculates deltaT
static double calculateTimeStep(const EulerSolverState *state){
    double maxSpeed = 0.0;
    for (unsigned int i = 0; i < state->config->nx; i++) {
        double rho = state->Q[i][0] / state->area[i];
        double u = state->Q[i][1] / state->Q[i][0];
        double p = pressure(state->Q[i], state->area[i], state->config->gamma);
        double a = sqrt(state->config->gamma * p / rho);
        maxSpeed = fmax(maxSpeed, fabs(u) + a);
    }
    return state->config->cfl * state->dx / maxSpeed;
}


void eulerInit(EulerSolverState *state, AlgorithmConfig *config){

    state->config = config;
    state->dx = (config->x1 - config->x0) / (config->nx - 1);
    state->time = 0.0;

    // sizeof(*state->Q) is the size of one row, so 3 doubles, not one

    state->x = calloc(config->nx, sizeof(*state->x));
    state->area = calloc(config->nx, sizeof(*state->area));
    state->darea_dx = calloc(config->nx, sizeof(*state->darea_dx));
    state->Q = calloc(config->nx, sizeof(*state->Q));
    state->Q_predictor = calloc(config->nx, sizeof(*state->Q_predictor));
    state->Q_next = calloc(config->nx, sizeof(*state->Q_next));

    //loop that fills the grid and the nozzle geometry
    for (unsigned int i = 0; i < config->nx; i++) {

        state->x[i] = config->x0 + (i * state->dx);

        // A(x):
        state->area[i] = 1.398 + 0.347 * tanh(0.8 * state->x[i] - 4.0);
        // A'(x):
        state->darea_dx[i] = 0.347 * 0.8 * (1.0 - tanh(0.8 * state->x[i] - 4.0) * tanh(0.8 * state->x[i] - 4.0));
    }

    eulerSetInitialCond(state);
    state->dt = calculateTimeStep(state);
}

void eulerClear(EulerSolverState *state){
    
    if (state == NULL) {
        return;
    }
    free(state->x);
    free(state->area);
    free(state->darea_dx);
    free(state->Q);
    free(state->Q_predictor);
    free(state->Q_next);
    *state = (EulerSolverState){0};
}

void eulerSetInitialCond(EulerSolverState *state){

    AlgorithmConfig *config = state->config;
    double rho = config->rho_in;
    double p = config->p_in;
    double u = config->mach_in * sqrt(config->gamma * p / rho);
    double e = p / (config->gamma - 1.0) + 0.5 * rho * u * u;

    // start with a uniform primitive state. This is only an initial guess;
    // the solution changes as we march towards the steady nozzle flow.
    for (unsigned int i = 0; i < config->nx; i++) {
        state->Q[i][0] = rho * state->area[i];
        state->Q[i][1] = rho * u * state->area[i];
        state->Q[i][2] = e * state->area[i];
    }

    setBoundaryCond(state, state->Q);

    // copy solution calculated above into the other Q arrays
    memcpy(state->Q_predictor, state->Q, config->nx * sizeof(*state->Q));
    memcpy(state->Q_next, state->Q, config->nx * sizeof(*state->Q));
}


static void calculateFlux(const EulerSolverState *state, double (*Q)[3], double (*F)[3], double (*S)[3])
{
    for (unsigned int i = 0; i < state->config->nx; i++) {

        double u = Q[i][1] / Q[i][0];
        double p = pressure(Q[i], state->area[i], state->config->gamma);

        // calculate F from given definition
        F[i][0] = Q[i][1];
        F[i][1] = Q[i][1] * u + p * state->area[i];
        F[i][2] = u * (Q[i][2] + p * state->area[i]);

        // calculate S from given definition
        S[i][0] = 0.0;
        S[i][1] = p * state->darea_dx[i];
        S[i][2] = 0.0;
    }
}

// the dissipation is added to the maccormack method to increase stability.
// TO DO: double check if this conflicts with anything written in the lecture notes
static void calculateDissipation(const EulerSolverState *state, double (*Q)[3], double (*dissipation)[3], double *coefficient)
{
    unsigned int n = state->config->nx;
    double *p = allocateArray(n, sizeof(*p));
    double *speed = allocateArray(n, sizeof(*speed));
    double *sensor = allocateArray(n, sizeof(*sensor));

    for (unsigned int i = 0; i < n; i++) {

        p[i] = pressure(Q[i], state->area[i], state->config->gamma);
        double rho = Q[i][0] / state->area[i];
        speed[i] = fabs(Q[i][1] / Q[i][0]) + sqrt(state->config->gamma * p[i] / rho);
    }

    for (unsigned int i = 1; i < n - 1; i++) {
        sensor[i] = fabs(p[i+1] - 2.0 * p[i] + p[i-1]) / (p[i+1] + 2.0 * p[i] + p[i-1]);
    }

    sensor[0] = sensor[1];
    sensor[n-1] = sensor[n-2];

    // dissipation[i] is a face flux between grid points i and i+1.
    // Using its difference in the update keeps the added term conservative.
    for (unsigned int i = 0; i < n - 1; i++) {

        double lambda = fmax(speed[i], speed[i+1]);
        double eps2 = DISSIPATION_2 * fmax(sensor[i], sensor[i+1]);
        double eps4 = fmax(0.0, DISSIPATION_4 - eps2);

        coefficient[i] = lambda * eps2;
        
        for (int k = 0; k < 3; k++) {
            double thirdDifference = 0.0;
            if (i > 0 && i < n - 2) {
                thirdDifference = Q[i+2][k] - 3.0 * Q[i+1][k] + 3.0 * Q[i][k] - Q[i-1][k];
            }
            dissipation[i][k] = lambda * (eps2 * (Q[i+1][k] - Q[i][k]) - eps4 * thirdDifference);
        }
    }
    free(p);
    free(speed);
    free(sensor);
}

// where the magic happens:
static void macCormackStep(EulerSolverState *state, double (*F)[3], double (*S)[3], double (*dissipation)[3], double *coefficient){
    unsigned int n = state->config->nx;
    double ratio = state->dt / state->dx;
    calculateFlux(state, state->Q, F, S);
    calculateDissipation(state, state->Q, dissipation, coefficient);

    // predictor: forward difference for the flux
    for (unsigned int i = 1; i < n - 1; i++) {
        for (int k = 0; k < 3; k++) {
            state->Q_predictor[i][k] =
            state->Q[i][k] - ratio * (F[i+1][k] - F[i][k]) + state->dt * S[i][k] + ratio * (dissipation[i][k] - dissipation[i-1][k]);
        }
    }

    // ensure boundary condition is met
    setBoundaryCond(state, state->Q_predictor);

    calculateFlux(state, state->Q_predictor, F, S);
    calculateDissipation(state, state->Q_predictor, dissipation, coefficient);

    // corrector: backward difference using the predicted state
    for (unsigned int i = 1; i < n - 1; i++) {
        for (int k = 0; k < 3; k++) {
            state->Q_next[i][k] = 0.5 * (state->Q[i][k] + state->Q_predictor[i][k]
                - ratio * (F[i][k] - F[i-1][k]) + state->dt * S[i][k]
                + ratio * (dissipation[i][k] - dissipation[i-1][k]));
        }
    }

    // ensure boundary condition is met, again
    setBoundaryCond(state, state->Q_next);
}

static void fluxJacobian(const double Q[3], double gamma, double J[3][3]){
    // J = dF/dQ. Area cancels here because Q and F both include it
    double u = Q[1] / Q[0];
    double H = gamma * Q[2] / Q[0] - 0.5 * (gamma - 1.0) * u * u;
    J[0][0] = 0.0;
    J[0][1] = 1.0;
    J[0][2] = 0.0;
    J[1][0] = 0.5 * (gamma - 3.0) * u * u;
    J[1][1] = (3.0 - gamma) * u;
    J[1][2] = gamma - 1.0;
    J[2][0] = u * (0.5 * (gamma - 1.0) * u * u - H);
    J[2][1] = H - (gamma - 1.0) * u * u;
    J[2][2] = gamma * u;
}


static void solveBlock(double B[3][3], double C[3][3], double rhs[3]){
    // solve B*[C_new,rhs_new] = [C,rhs] using Gaussian elimination.
    // Partial pivoting avoids dividing by a small diagonal when another row
    // provides a better pivot. B is overwritten with the identity matrix.
    for (int k = 0; k < 3; k++) {
        int pivot = k;
        for (int r = k + 1; r < 3; r++) {
            if (fabs(B[r][k]) > fabs(B[pivot][k])) {
                pivot = r;
            }
        }
        
        if (pivot != k) {
            for (int c = 0; c < 3; c++) {
                double temp = B[k][c]; B[k][c] = B[pivot][c]; B[pivot][c] = temp;
                temp = C[k][c]; C[k][c] = C[pivot][c]; C[pivot][c] = temp;
            }
            double temp = rhs[k]; rhs[k] = rhs[pivot]; rhs[pivot] = temp;
        }
        double diagonal = B[k][k];
        for (int c = 0; c < 3; c++) {
            B[k][c] /= diagonal;
            C[k][c] /= diagonal;
        }
        rhs[k] /= diagonal;
        for (int r = 0; r < 3; r++) {
            if (r == k) {
                continue;
            }
            double factor = B[r][k];
            for (int c = 0; c < 3; c++) {
                B[r][c] -= factor * B[k][c];
                C[r][c] -= factor * C[k][c];
            }
            rhs[r] -= factor * rhs[k];
        }
    }
}

// see if I cannot just reuse the function from the wavesolver
static void solveBlockThomas(double (*lower)[3][3], double (*diagonal)[3][3], double (*upper)[3][3], double (*rhs)[3], unsigned int n){
    // same idea as the scalar Thomas algorithm in Q1, but each entry is
    // now a 3x3 matrix because mass, momentum and energy are coupled.
    solveBlock(diagonal[0], upper[0], rhs[0]);
    for (unsigned int i = 1; i < n; i++) {
        for (int r = 0; r < 3; r++) {
            for (int k = 0; k < 3; k++) {
                rhs[i][r] -= lower[i][r][k] * rhs[i-1][k];
                for (int c = 0; c < 3; c++) {
                    diagonal[i][r][c] -= lower[i][r][k] * upper[i-1][k][c];
                }
            }
        }
        solveBlock(diagonal[i], upper[i], rhs[i]);
    }
    // back substitution. rhs gets overwritten with the solution delta Q.
    for (unsigned int i = n - 1; i > 0; i--) {
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) {
                rhs[i-1][r] -= upper[i-1][r][c] * rhs[i][c];
            }
        }
    }
}

static int beamWarmingStep(EulerSolverState *state, double (*F)[3], double (*S)[3], double (*dissipation)[3], double *coefficient){
    unsigned int n = state->config->nx;
    double gamma = state->config->gamma;
    double ratio = state->dt / state->dx;
    double (*lower)[3][3] = allocateArray(n, sizeof(*lower));
    double (*diagonal)[3][3] = allocateArray(n, sizeof(*diagonal));
    double (*upper)[3][3] = allocateArray(n, sizeof(*upper));
    double (*J)[3][3] = allocateArray(n, sizeof(*J));
    double (*delta)[3] = state->Q_predictor;

    calculateFlux(state, state->Q, F, S);
    calculateDissipation(state, state->Q, dissipation, coefficient);
    memset(delta, 0, n * sizeof(*delta));
    for (unsigned int i = 0; i < n; i++) {
        fluxJacobian(state->Q[i], gamma, J[i]);
        for (int k = 0; k < 3; k++) {
            diagonal[i][k][k] = 1.0;
        }
    }

    // Beam-Warming delta form, using backward Euler in time (theta = 1):
    // [I + dt*Dx*J - dt*dS/dQ - dt*D2] deltaQ = dt*(-Dx*F + S + D).
    // Dx is centered. Second-difference dissipation is implicit with frozen
    // coefficients; the fourth-difference part is explicit on the RHS.
    // This is first order in time, appropriate here for marching to steady state.
    for (unsigned int i = 1; i < n - 1; i++) {
        for (int r = 0; r < 3; r++) {
            delta[i][r] = -0.5 * ratio * (F[i+1][r] - F[i-1][r])
                + state->dt * S[i][r]
                + ratio * (dissipation[i][r] - dissipation[i-1][r]);
            for (int c = 0; c < 3; c++) {
                lower[i][r][c] = -0.5 * ratio * J[i-1][r][c];
                upper[i][r][c] = 0.5 * ratio * J[i+1][r][c];
            }
            lower[i][r][r] -= ratio * coefficient[i-1];
            diagonal[i][r][r] += ratio * (coefficient[i-1] + coefficient[i]);
            upper[i][r][r] -= ratio * coefficient[i];
        }
        // only the momentum source is nonzero: S[1] = p*dA/dx
        double u = state->Q[i][1] / state->Q[i][0];
        double factor = state->dt * state->darea_dx[i] * (gamma - 1.0) / state->area[i];
        diagonal[i][1][0] -= factor * 0.5 * u * u;
        diagonal[i][1][1] += factor * u;
        diagonal[i][1][2] -= factor;
    }

    // first row: deltaQ = 0 at the fixed inlet.
    // last row: linearize the outlet extrapolation, including fixed p_out
    // for the subsonic case, instead of freezing all outlet quantities.
    double areaRatio = state->area[n-1] / state->area[n-2];
    for (int k = 0; k < 3; k++) {
        lower[n-1][k][k] = -areaRatio;
    }
    if (state->config->outlet_type == SUBSONIC_OUTLET) {
        double u = state->Q[n-2][1] / state->Q[n-2][0];
        lower[n-1][2][0] = 0.5 * areaRatio * u * u;
        lower[n-1][2][1] = -areaRatio * u;
        lower[n-1][2][2] = 0.0;
    }

    solveBlockThomas(lower, diagonal, upper, delta, n);
    for (unsigned int i = 0; i < n; i++) {
        for (int k = 0; k < 3; k++) {
            state->Q_next[i][k] = state->Q[i][k] + delta[i][k];
        }
    }
    setBoundaryCond(state, state->Q_next);
    free(lower);
    free(diagonal);
    free(upper);
    free(J);
    return validSolution(state, state->Q_next);
}

// def some simplification necessary here
void eulerStep(EulerSolverState *state){
    unsigned int n = state->config->nx;
    setBoundaryCond(state, state->Q);
    
    state->dt = calculateTimeStep(state);
    double (*F)[3] = allocateArray(n, sizeof(*F));
    double (*S)[3] = allocateArray(n, sizeof(*S));
    double (*dissipation)[3] = allocateArray(n, sizeof(*dissipation));
    double *coefficient = allocateArray(n, sizeof(*coefficient));

    
    int accepted = 0;
    
    switch (state->config->algorithm) {
        case MACCORMACK:
            macCormackStep(state, F, S, dissipation, coefficient);
            break;

        case BEAM_WARMING:
            accepted = beamWarmingStep(state, F, S, dissipation, coefficient);
            break;

        default:
            break;
    }
    
    free(F);
    free(S);
    free(dissipation);
    free(coefficient);

    // increment time and swap the solution arrays, same as in Q1
    state->time += state->dt;
    double (*temp)[3] = state->Q;
    state->Q = state->Q_next;
    state->Q_next = temp;
}


void stateWriteToCSV(FILE *outputFile, const EulerSolverState *state){

    if (outputFile == NULL) {
        solverError("CSV file could not be opened");
    }

    // one row per grid point. Calling this again appends another snapshot.
    // Open the file with "w" for a new run; write the header only once.
    if (ftell(outputFile) == 0) {
        fprintf(outputFile, "t,x,A,rho,u,p,M,e,mass_flow\n");
    }
    
    for (unsigned int i = 0; i < state->config->nx; i++) {
        double rho = state->Q[i][0] / state->area[i];
        double u = state->Q[i][1] / state->Q[i][0];
        double p = pressure(state->Q[i], state->area[i], state->config->gamma);
        double a = sqrt(state->config->gamma * p / rho);
        double e = state->Q[i][2] / state->area[i];
        fprintf(outputFile, "%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
                state->time, state->x[i], state->area[i], rho, u, p, u / a, e, state->Q[i][1]);
    }
    if (ferror(outputFile)) {
        solverError("could not write CSV output");
    }
}

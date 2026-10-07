#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "eulersolver.h"



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
        unsigned int inside = last - 1;
        double gamma = config->gamma;

        // get the flow conditions at the last interior point
        double rho_inside = Q[inside][0] / state->area[inside];
        double u_inside = Q[inside][1] / Q[inside][0];
        double p_inside = pressure(Q[inside], state->area[inside], gamma);

        // prescribed outlet pressure
        double p_out = config->back_pressure_ratio * config->p_in;

        // preserve the entropy from the interior: p / rho^gamma = constant
        double rho_out = rho_inside * pow(p_out / p_inside, 1.0 / gamma);

        // speed of sound at the interior point and at the outlet
        double c_inside = sqrt(gamma * p_inside / rho_inside);
        double c_out = sqrt(gamma * p_out / rho_out);

        // preserve the right-running invariant: u + 2c/(gamma-1)
        double u_out = u_inside
                    + 2.0 * (c_inside - c_out) / (gamma - 1.0);

        // convert the outlet flow conditions back to conserved quantities
        double e_out = p_out / (gamma - 1.0)
                    + 0.5 * rho_out * u_out * u_out;

        Q[last][0] = rho_out * state->area[last];
        Q[last][1] = rho_out * u_out * state->area[last];
        Q[last][2] = e_out * state->area[last];
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


// where the magic happens:
static void macCormackStep(EulerSolverState *state, double (*F)[3], double (*S)[3]){
    unsigned int n = state->config->nx;
    double ratio = state->dt / state->dx;
    calculateFlux(state, state->Q, F, S);


    // predictor: forward difference for the flux
    for (unsigned int i = 1; i < n - 1; i++) {
        for (int k = 0; k < 3; k++) {
            state->Q_predictor[i][k] =
            state->Q[i][k] - ratio * (F[i+1][k] - F[i][k]) + state->dt * S[i][k];
        }
    }

    // ensure boundary condition is met
    setBoundaryCond(state, state->Q_predictor);

    calculateFlux(state, state->Q_predictor, F, S);

    // corrector: backward difference using the predicted state
    for (unsigned int i = 1; i < n - 1; i++) {
        for (int k = 0; k < 3; k++) {
            state->Q_next[i][k] = 
            0.5 * (state->Q[i][k] + state->Q_predictor[i][k] - ratio * (F[i][k] - F[i-1][k]) + state->dt * S[i][k]);
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

static void beamWarmingStep(EulerSolverState *state, double (*F)[3], double (*S)[3]){
    
    unsigned int n = state->config->nx;
    double gamma = state->config->gamma;
    double ratio = state->dt / state->dx;
    
    double (*lower)[3][3] = calloc(n, sizeof(*lower));
    double (*diagonal)[3][3] = calloc(n, sizeof(*diagonal));
    double (*upper)[3][3] = calloc(n, sizeof(*upper));
    double (*J)[3][3] = calloc(n, sizeof(*J));
    double (*delta)[3] = state->Q_predictor;

    double dissipation_explicit = 0.05;
    double dissipation_implicit = 0.15;

    double eps_e = state->dt * dissipation_explicit;
    double eps_i = state->dt * dissipation_implicit;

    calculateFlux(state, state->Q, F, S);

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

            delta[i][r] = -0.5 * ratio * (F[i+1][r] - F[i-1][r]) + state->dt * S[i][r];

            for (int c = 0; c < 3; c++) {
                lower[i][r][c] = -0.5 * ratio * J[i-1][r][c];
                upper[i][r][c] = 0.5 * ratio * J[i+1][r][c];
            }

            // implicit second difference acting on delta Q
            lower[i][r][r] -= eps_i;
            diagonal[i][r][r] += 2.0 * eps_i;
            upper[i][r][r] -= eps_i;

            // explicit fourth difference acting on the current solution.
            // The full stencil is available only at points 2 through n-3.
            if (i >= 2 && i < n - 2) {
                double fourth_difference =
                    state->Q[i-2][r]
                    - 4.0 * state->Q[i-1][r]
                    + 6.0 * state->Q[i][r]
                    - 4.0 * state->Q[i+1][r]
                    +       state->Q[i+2][r];

                delta[i][r] -= eps_e * fourth_difference;
            }
        }

        // only the momentum source is nonzero: S[1] = p*dA/dx
        double u = state->Q[i][1] / state->Q[i][0];
        double factor = state->dt * state->darea_dx[i] * (gamma - 1.0) / state->area[i];
        diagonal[i][1][0] -= factor * 0.5 * u * u;
        diagonal[i][1][1] += factor * u;
        diagonal[i][1][2] -= factor;
    }

    
    unsigned int last = n - 1;
    unsigned int inside = n - 2;

    double A_inside = state->area[inside];
    double A_out = state->area[last];

    if (state->config->outlet_type == SUPERSONIC_OUTLET) {
        // same extrapolation as before
        double areaRatio = A_out / A_inside;

        for (int k = 0; k < 3; k++) {
            lower[last][k][k] = -areaRatio;
        }
    } else {
        // recover the interior state
        double rho = state->Q[inside][0] / A_inside;
        double u = state->Q[inside][1] / state->Q[inside][0];
        double p = pressure(state->Q[inside], A_inside, gamma);

        // calculate the same outlet state as setBoundaryCond()
        double p_out = state->config->back_pressure_ratio
                    * state->config->p_in;

        double densityRatio = pow(p_out / p, 1.0 / gamma);
        double rho_out = rho * densityRatio;

        double c = sqrt(gamma * p / rho);
        double c_out = sqrt(gamma * p_out / rho_out);
        double u_out = u + 2.0 * (c - c_out) / (gamma - 1.0);

        // derivatives of interior rho, u and p with respect to Q[0], Q[1], Q[2]
        double drho[3] = {
            1.0 / A_inside, 0.0, 0.0
        };

        double du[3] = {
            -u / state->Q[inside][0],
            1.0 / state->Q[inside][0],
            0.0
        };

        double dp[3] = {
            (gamma - 1.0) * 0.5 * u * u / A_inside,
            -(gamma - 1.0) * u / A_inside,
            (gamma - 1.0) / A_inside
        };

        // build each column of the boundary Jacobian using the chain rule
        // p_out is prescribed, so its derivative is zero
        for (int k = 0; k < 3; k++) {
            double drho_out = densityRatio * drho[k]
                        - rho_out * dp[k] / (gamma * p);

            double dc = 0.5 * c * (dp[k] / p - drho[k] / rho);
            double dc_out = -0.5 * c_out * drho_out / rho_out;

            double du_out = du[k]
                        + 2.0 * (dc - dc_out) / (gamma - 1.0);

            // lower = -dQ_out/dQ_inside
            lower[last][0][k] = -A_out * drho_out;

            lower[last][1][k] = -A_out
                            * (u_out * drho_out + rho_out * du_out);

            lower[last][2][k] = -A_out
                            * (0.5 * u_out * u_out * drho_out
                                + rho_out * u_out * du_out);
        }
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

}

// def some simplification necessary here
void eulerStep(EulerSolverState *state){
    unsigned int n = state->config->nx;
    setBoundaryCond(state, state->Q);
    
    state->dt = calculateTimeStep(state);
    double (*F)[3] = calloc(n, sizeof(*F));
    double (*S)[3] = calloc(n, sizeof(*S));

    
    switch (state->config->algorithm) {
        case MACCORMACK:
            macCormackStep(state, F, S);
            break;

        case BEAM_WARMING:
            beamWarmingStep(state, F, S);
            break;

        default:
            break;
    }
    
    free(F);
    free(S);

    // increment time and swap the solution arrays, same as in Q1
    state->time += state->dt;
    double (*temp)[3] = state->Q;
    state->Q = state->Q_next;
    state->Q_next = temp;
}


void stateWriteToCSV(FILE *outputFile, const EulerSolverState *state){

    if (outputFile == NULL) {
        printf("CSV file could not be opened");
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
        printf("could not write CSV output");
    }
}

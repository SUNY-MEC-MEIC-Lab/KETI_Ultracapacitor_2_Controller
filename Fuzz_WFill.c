// This is a project to generate the waterfill function with inputs of SoA, k_nominal, and dt

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Fuzz_WFill.h"

//int main() { //EXAMPLE OF FUZZY CONTROLLER CODE WITH 4 CELLS, SHOULD BE USED WITH 48
//    int N = 5;
//    double SoA[5] = { 0.8, 0.85, 0.75, 0.9, 0.82 };
//    double aging_rate[5] = { 0.02, 0.025, 0.03, 0.015, 0.02 };
//    double F_gain = 0.5;
//    double u[5];
//
//    fuzzy_controller(SoA, N, aging_rate, F_gain, u);
//
//    for (int i = 0; i < N; i++) {
//        printf("u[%d] = %f\n", i, u[i]);
//    }
//
//    return 0;
//}

//int main() { //EXAMPLE OF LINEAR WATERFILL CONTROLLER CODE WITH 48. SETUP FOR USE WITH 48 CELLS ONLY
//    int N = 24;
//    double SoA[] = { 0.8, 0.7, 0.6, 0.9, ... };
//    double k_nominal[] = { 0.01, 0.02, 0.015, 0.017, ... };
//    double dt = 1.0;
//    double u[24];
//
//    linear_waterfill_controller(N, SoA, k_nominal, dt, u);
//
//    for (int i = 0; i < N; i++) {
//        printf("u[%d] = %f\n", i, u[i]);
//    }
//
//    return 0;
//}


void fuzzy_controller(double* SoA, int N, double* aging_rate, double F_gain, double* u) {
    double soa_sum = 0.0, aging_rate_sum = 0.0;
    for (int i = 0; i < N; i++) {
        soa_sum += SoA[i];
        aging_rate_sum += aging_rate[i];
    }

    double soa_avg = soa_sum / N;
    double aging_rate_avg = aging_rate_sum / N;

    // Compute errors and their standard deviations
    double* soa_error = (double*)malloc(N * sizeof(double));
    double* aging_rate_error = (double*)malloc(N * sizeof(double));
    double soa_var = 0.0, rate_var = 0.0;

    for (int i = 0; i < N; i++) {
        soa_error[i] = SoA[i] - soa_avg;
        aging_rate_error[i] = aging_rate[i] - aging_rate_avg;
        soa_var += soa_error[i] * soa_error[i];
        rate_var += aging_rate_error[i] * aging_rate_error[i];
    }

    double soa_thresh = sqrt(soa_var / N);
    double rate_thresh = sqrt(rate_var / N);
    if (rate_thresh < 1e-9) rate_thresh = 1e-9;

    // Fuzzy adjustment
    int* adj = (int*)calloc(N, sizeof(int));  // initialize with zeros

    for (int i = 0; i < N; i++) {
        if (soa_error[i] > soa_thresh) adj[i] += 1;
        if (soa_error[i] < -soa_thresh) adj[i] -= 1;
        if (aging_rate_error[i] > rate_thresh) adj[i] += 1;
        if (aging_rate_error[i] < -rate_thresh) adj[i] -= 1;
    }

    // Compute u[i] = (1/N) - F_gain * adj[i]/N
    double sum_u = 0.0;
    for (int i = 0; i < N; i++) {
        u[i] = (1.0 / N) - (F_gain * adj[i] / N);
        if (u[i] < 0.0) u[i] = 0.0;
        sum_u += u[i];
    }

    // Normalize u
    if (sum_u > 0.0) {
        for (int i = 0; i < N; i++) {
            u[i] /= sum_u;
        }
    }

    // Cleanup
    free(soa_error);
    free(aging_rate_error);
    free(adj);
}

void linear_waterfill_controller(int nc, double SoA [], double k_nominal [], double dt, double* u) { // Waterfill controller (linear model)

    double lo = 0.0;
    double hi = 0.0;
    double mid = 0.0;
    double s = 0.0;
    double u_tot = 0.0;

    lo = my_max(SoA, nc);
    hi = lo + mean(k_nominal, nc) * dt;

    while (total_load(hi, SoA, k_nominal, dt) < 1) { hi = hi + (hi - lo); }

    for (int i = 1; i <= 1; i = i + 1) {
        mid = (lo + hi) / 2.0;
        if (total_load(mid, SoA, k_nominal, dt) > 1.0)
            hi = mid;
        else
            lo = mid;
    }
    s = (lo + hi) / 2.0;
    for (int j = 0; j < 48; j = j + 1) {
        u[j] = (s - SoA[j]) / (k_nominal[j] * dt);
    }
    for (int k = 0; k < 48; k++) {
        u[k] = my_min_el(1, u[k]);
    }
    for (int l = 0; l < 48; l++) {
        u[l] = my_max_el(0, u[l]);
    }
    for (int m = 0; m < 48; m++) {
        u_tot = u_tot + u[m];
    }
    for (int n = 0; n < 48; n++) {
        u[n] = u[n] / u_tot;
    }
    // return *u;
}

double total_load(double target_s, double SoA [], double k_nominal [], double dt) {
    double u_trial[48];
    double u_clamped[48];
    double min_uc[48];
    double g = 0.0;
    for (int i = 0; i < 48; i++) {
        u_trial[i] = (target_s - SoA[i]) / (k_nominal[i] * dt);
    }
    for (int i = 0; i < 48; i++) {
        min_uc[i] = my_min_el(1.0,u_trial[i]);
    }
    for (int i = 0; i < 48; i++) {
        u_clamped[i] = my_max_el(0, min_uc[i]);
    }
    for (int i = 0; i < 48; i++) {
        g = g + u_clamped[i];
    }
    return g;
}

double my_max(double vectorset [], int vecsize) {
    double max_val = vectorset[0];
    for (int i = 0; i < vecsize; i = i + 1) {
        if (max_val < vectorset[i])
            max_val = vectorset[i];
    }
    return max_val;
}

double my_max_el(double ele, double vec_el) {
    if (ele > vec_el)
        vec_el = ele;
    return vec_el;
}

double my_min_el(double ele, double vec_el) {
    if (ele < vec_el)
        vec_el = ele;
    return vec_el;
}

double my_min(double vectorset[], int vecsize) {
    double min_val = vectorset[0];
    for (int i = 0; i < vecsize; i = i + 1) {
        if (min_val > vectorset[i])
            min_val = vectorset[i];
    }
    return min_val;
}

double mean(double vectorset[], int vecsize) {
    double mean_val = 0;
    for (int i = 0; i < vecsize; i = i + 1) {
        mean_val = mean_val + vectorset[i];
    }
    mean_val = mean_val / vecsize;
    return mean_val;
}
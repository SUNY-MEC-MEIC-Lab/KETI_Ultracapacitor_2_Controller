#pragma once
#ifndef Fuzz_WFill_contlr_H

#define Fuzz_WFill_contlr_H

double total_load(double target_s, double SoA[], double k_nominal[], double dt);

void linear_waterfill_controller(const int nc, double SoA[], double k_nominal[], double dt, double* u);

void fuzzy_controller(double* SoA, int N, double* aging_rate, double F_gain, double* u);

double my_max(double vectorset[], int vecsize);

double my_max_el(double ele, double vec_el);

double my_min_el(double ele, double vec_el);

double my_min(double vectorset[], int vecsize);

double mean(double vectorset[], int vecsize);

#endif // Fuzz_WFill_contlr_H
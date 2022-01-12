#ifndef gshape_h
#define gshape_h
#include <fstream>
#include <iostream>

//
double LSolver6(double hess[6][6], double der6[6], double dx[6]);

double cal_cart_dervs(int na1, double* xyz1, double* w1, int na2, double* xyz2,
                         double* w2, double* cderv, double* chess);
double cal_cart_dervs2(int na1, double* xyz1, double* alpha1, double* w1, int na2, double* xyz2,
                         double* alpha2, double* w2, double* cderv, double* chess, double* wxpp, 
                         double* afac, double* ahaf);

void cart_dervs_2rt(int na1, double* xyz1, double* cderv, double* chess, double rt_derv[6],
                       double rt_hess[6][6]);

#endif /* gshape_h */

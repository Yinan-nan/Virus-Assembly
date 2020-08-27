
#ifndef HIV2_PARAMETERS_H
#define HIV2_PARAMETERS_H


 //HIV2_PARAMETERS_H

/*class parameters {
public:
    double kb;
    double ks;
    double th0;
    double sigma = .0001;
*//*double r0;*//*
#define R3_2 0.86602540378443864676
    int MAX_num_tri = 200;
};*/

// #define ks sqrt(0.2)
// #define kb (1/ks)
#define ks sqrt(18)
#define kb (1/ks)
#define sigma  .0001
// #define r0 1.12
//#define r0 1.32
#define r0 1.6
#define R3_2 0.86602540378443864676
#define MAX_num_tri 600
#define max_iter 300
#define epsi 1.0
// #define sigma 1.0
#define pi 3.1415926
#define ri 0.951057
#define b0 1


#endif


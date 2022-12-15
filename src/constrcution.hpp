#include <cmath>
#include <iostream>
#include "define.hpp"

class CSSWM {
public:
    class patch {
    public: 
        // constructer
        patch();

        double hp[NX][NY], h[NX][NY], hm[NX][NY];
        double up[NX][NY], u[NX][NY], um[NX][NY];
        double vp[NX][NY], v[NX][NY], vm[NX][NY];

        double lon[NX][NY], lat[NX][NY];

        double x[NX][NY], y[NX][NY];

        double A[NX][NY][4], IA[NX][NY][4];
    };

    CSSWM();
    patch csswm[6];
    double sqrtG[NX][NY], gamma[NX][NY], gLower[NX][NY][4], gUpper[NX][NY][4];

    // ***********************************************************************************
    // In transform.cpp
    double Cube2Sphere_U(CSSWM &, int, int, int);
    double Cube2Sphere_V(CSSWM &, int, int, int);
    double Sphere2Cube_U(CSSWM &, int, int, int);
    double Sphere2Cube_V(CSSWM &, int, int, int);
    double Cube2Cube_U(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    double Cube2Cube_V(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    double Cube2Cube_BV2AU(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    double Cube2Cube_BU2AV(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    void matrixMul(double firstMatrix[4], double secondMatrix[4], double mult[2][2]);
    // ***********************************************************************************

    // ***********************************************************************************
    // In bp_h.cpp
    void BP_h(CSSWM &);
    double interpolate(double A1, double A2, double V1, double V2, double B);
    // ***********************************************************************************


private:
    void Construct_gamma_sqrtG_GUpper(double **alpha2D, double **beta2D, double gamma[NX][NY], double sqrtG[NX][NY], double gUpper[NX][NY][4], double gLower[NX][NY][4]);
    void Construct_p0123_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]);
    void Construct_p4_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]);
    void Construct_p5_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]);
};
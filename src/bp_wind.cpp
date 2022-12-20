#include "constrcution.hpp"

void CSSWM::BP_wind_convert(CSSWM &model) {
    for (int idx = 0; idx < NX; idx++) {
        // patch0 (Left, Right, Up, Down)
        model.csswm[0].up[0][idx] = model.Cube2Cube_U(model, 0, 3, 0, idx, NX-2, idx);
        model.csswm[0].vp[0][idx] = model.Cube2Cube_V(model, 0, 3, 0, idx, NX-2, idx);

        model.csswm[0].up[NX-1][idx] = model.Cube2Cube_U(model, 0, 1, NX-1, idx, 1, idx);
        model.csswm[0].vp[NX-1][idx] = model.Cube2Cube_V(model, 0, 1, NX-1, idx, 1, idx);

        model.csswm[0].up[idx][NY-1] = model.Cube2Cube_U(model, 0, 4, idx, NY-1, idx, 1);
        model.csswm[0].vp[idx][NY-1] = model.Cube2Cube_V(model, 0, 4, idx, NY-1, idx, 1);

        model.csswm[0].up[idx][0] = model.Cube2Cube_U(model, 0, 5, idx, 0, idx, NY-2);
        model.csswm[0].vp[idx][0] = model.Cube2Cube_V(model, 0, 5, idx, 0, idx, NY-2);

        // patch1 (Left, Right, Up, Down)
        model.csswm[1].up[0][idx] = model.Cube2Cube_U(model, 1, 0, 0, idx, NX-2, idx);
        model.csswm[1].vp[0][idx] = model.Cube2Cube_V(model, 1, 0, 0, idx, NX-2, idx);

        model.csswm[1].up[NX-1][idx] = model.Cube2Cube_U(model, 1, 2, NX-1, idx, 1, idx);
        model.csswm[1].vp[NX-1][idx] = model.Cube2Cube_V(model, 1, 2, NX-1, idx, 1, idx);

        model.csswm[1].up[idx][NY-1] = model.Cube2Cube_BV2AU(model, 1, 4, idx, NY-1, NX-2, idx);
        model.csswm[1].vp[idx][NY-1] = model.Cube2Cube_BU2AV(model, 1, 4, idx, NY-1, NX-2, idx);

        model.csswm[1].up[idx][0] = model.Cube2Cube_BV2AU(model, 1, 5, idx, 0, NX-2, NY-1-idx);
        model.csswm[1].vp[idx][0] = model.Cube2Cube_BU2AV(model, 1, 5, idx, 0, NX-2, NY-1-idx);

        // patch2 (Left, Right, Up, Down)
        model.csswm[2].up[0][idx] = model.Cube2Cube_U(model, 2, 1, 0, idx, NX-2, idx);
        model.csswm[2].vp[0][idx] = model.Cube2Cube_V(model, 2, 1, 0, idx, NX-2, idx);

        model.csswm[2].up[NX-1][idx] = model.Cube2Cube_U(model, 2, 3, NX-1, idx, 1, idx);
        model.csswm[2].vp[NX-1][idx] = model.Cube2Cube_V(model, 2, 3, NX-1, idx, 1, idx);

        model.csswm[2].up[idx][NY-1] = model.Cube2Cube_U(model, 2, 4, idx, NY-1, NX-1-idx, NY-2);
        model.csswm[2].vp[idx][NY-1] = model.Cube2Cube_V(model, 2, 4, idx, NY-1, NX-1-idx, NY-2);

        model.csswm[2].up[idx][0] = model.Cube2Cube_U(model, 2, 5, idx, 0, NX-1-idx, 1);
        model.csswm[2].vp[idx][0] = model.Cube2Cube_V(model, 2, 5, idx, 0, NX-1-idx, 1);

        // patch3 (Left, Right, Up, Down)
        model.csswm[3].up[0][idx] = model.Cube2Cube_U(model, 3, 2, 0, idx, NX-2, idx);
        model.csswm[3].vp[0][idx] = model.Cube2Cube_V(model, 3, 2, 0, idx, NX-2, idx);

        model.csswm[3].up[NX-1][idx] = model.Cube2Cube_U(model, 3, 0, NX-1, idx, 1, idx);
        model.csswm[3].vp[NX-1][idx] = model.Cube2Cube_V(model, 3, 0, NX-1, idx, 1, idx);

        model.csswm[3].up[idx][NY-1] = model.Cube2Cube_BV2AU(model, 3, 4, idx, NY-1, 1, NY-1-idx);
        model.csswm[3].vp[idx][NY-1] = model.Cube2Cube_BU2AV(model, 3, 4, idx, NY-1, 1, NY-1-idx);

        model.csswm[3].up[idx][0] = model.Cube2Cube_BV2AU(model, 3, 5, idx, 0, 1, idx);
        model.csswm[3].vp[idx][0] = model.Cube2Cube_BU2AV(model, 3, 5, idx, 0, 1, idx);

        // patch4 (Left, Right, Up, Down)
        model.csswm[4].up[0][idx] = model.Cube2Cube_BV2AU(model, 4, 3, 0, idx, NX-1-idx, NY-2);
        model.csswm[4].vp[0][idx] = model.Cube2Cube_BU2AV(model, 4, 3, 0, idx, NX-1-idx, NY-2);

        model.csswm[4].up[NX-1][idx] = model.Cube2Cube_BV2AU(model, 4, 1, NX-1, idx, idx, NY-2);
        model.csswm[4].vp[NX-1][idx] = model.Cube2Cube_BU2AV(model, 4, 1, NX-1, idx, idx, NY-2);

        model.csswm[4].up[idx][NY-1] = model.Cube2Cube_U(model, 4, 2, idx, NY-1, NX-1-idx, NY-2);
        model.csswm[4].vp[idx][NY-1] = model.Cube2Cube_V(model, 4, 2, idx, NY-1, NX-1-idx, NY-2);

        model.csswm[4].up[idx][0] = model.Cube2Cube_U(model, 4, 0, idx, 0, idx, NY-2);
        model.csswm[4].vp[idx][0] = model.Cube2Cube_V(model, 4, 0, idx, 0, idx, NY-2);

        // patch5 (Left, Right, Up, Down)
        model.csswm[5].up[0][idx] = model.Cube2Cube_BV2AU(model, 5, 3, 0, idx, idx, 1);
        model.csswm[5].vp[0][idx] = model.Cube2Cube_BU2AV(model, 5, 3, 0, idx, idx, 1);

        model.csswm[5].up[NX-1][idx] = model.Cube2Cube_BV2AU(model, 5, 1, NX-1, idx, NX-1-idx, 1);
        model.csswm[5].vp[NX-1][idx] = model.Cube2Cube_BU2AV(model, 5, 1, NX-1, idx, NX-1-idx, 1);

        model.csswm[5].up[idx][NY-1] = model.Cube2Cube_U(model, 5, 0, idx, NY-1, idx, 1);
        model.csswm[5].vp[idx][NY-1] = model.Cube2Cube_V(model, 5, 0, idx, NY-1, idx, 1);

        model.csswm[5].up[idx][0] = model.Cube2Cube_U(model, 5, 2, idx, 0, NX-1-idx, 1);
        model.csswm[5].vp[idx][0] = model.Cube2Cube_V(model, 5, 2, idx, 0, NX-1-idx, 1);
    }
}


void CSSWM::BP_wind_interpolation(CSSWM &model) {
    double B, A1, A2, V1, V2, V3, V4;
    int idx = 0, ipIdx = 0;
    int checkIP[NX][2];

    // Construct Interpolation 2D array filled with index (Note that all interpolation relation between every boundary is same)
    while (idx < NX && ipIdx < NX - 1) {
        B = model.csswm[0].lat[NX-1][idx];
        A1 = model.csswm[1].lat[1][ipIdx], A2 = model.csswm[1].lat[1][ipIdx+1];

        if (A1 < B && B < A2) {
            checkIP[idx][0] = ipIdx;
            checkIP[idx][1] = ipIdx + 1;
            idx++;
        }
        else if (A1 == B) {
            checkIP[idx][0] = ipIdx;
            checkIP[idx][1] = ipIdx;
            idx++;
        }
        else if (A2 == B) {
            checkIP[idx][0] = ipIdx + 1;
            checkIP[idx][1] = ipIdx + 1;
            idx++;
        }
        else {
            ipIdx++;
        }
    }

    // Construct a array for dealing with interpolation between all patch
    // (p1, p2, i1, j1, i2, j2, reversed, LonLat: 1/0) 
    // Left, Right, Up, Down
    int match[24][8] = {
        {0, 3, 0, -1, NX-2, -1, 0, 0},  {0, 1, NX-1, -1, 1, -1, 0, 0},    {0, 4, -1, NY-1, -1, 1, 0, 1},     {0, 5, -1, 0, -1, NY-2, 0, 1},
        {1, 0, 0, -1, NX-2, -1, 0, 0},  {1, 2, NX-1, -1, 1, -1, 0, 0},    {1, 4, -1, NY-1, NX-2, -1, 0, 1},  {1, 5, -1, 0, NX-2, -1, 1, 1},
        {2, 1, 0, -1, NX-2, -1, 0, 0},  {2, 3, NX-1, -1, 1, -1, 0, 0},    {2, 4, -1, NY-1, -1, NY-2, 1, 1},  {2, 5, -1, 0, -1, 1, 1, 1},
        {3, 2, 0, -1, NX-2, -1, 0, 0},  {3, 0, NX-1, -1, 1, -1, 0, 0},    {3, 4, -1, NY-1, 1, -1, 1, 1},     {3, 5, -1, 0, 1, -1, 0, 1},
        {4, 3, 0, -1, -1, NY-2, 1, 1},  {4, 1, NX-1, -1, -1, NY-2, 0, 1}, {4, 2, -1, NY-1, -1, NY-2, 1, 1},  {4, 0, -1, 0, -1, NY-2, 0, 1},
        {5, 3, 0, -1, -1, 1, 0, 1},     {5, 1, NX-1, -1, -1, 1, 1, 1},    {5, 0, -1, NY-1, -1, 1, 0, 1},     {5, 2, -1, 0, -1, 1, 1, 1}
    };

    int p1, p2, i1, j1, i2, j2, reversed, lonlat;
    double uIP, vIP;
    double alpha, beta;
    double alpha_B, beta_B;
    double gLower[4], IA[4], A[4], gUpper[4];
    for (int pp = 0; pp < 24; pp++) {
        p1 = match[pp][0], p2 = match[pp][1], i1 = match[pp][2], j1 = match[pp][3], i2 = match[pp][4], j2 = match[pp][5], reversed = match[pp][6], lonlat = match[pp][7];
        for (int idx = 0; idx < NX; idx++) {
            if (lonlat == 0) {
                int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                int I2_1 = i2 == -1 ? reversed ? checkIP[NX-1-idx][0] : checkIP[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? checkIP[NY-1-idx][0] : checkIP[idx][0] : j2;
                int I2_2 = i2 == -1 ? reversed ? checkIP[NX-1-idx][1] : checkIP[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? checkIP[NY-1-idx][1] : checkIP[idx][1] : j2;

                B = model.csswm[p1].lat[I1][J1];
                A1 = model.csswm[p2].lat[I2_1][J2_1], A2 = model.csswm[p2].lat[I2_2][J2_2];
                V1 = model.csswm[p2].up[I2_1][J2_1], V2 = model.csswm[p2].up[I2_2][J2_2];
                V3 = model.csswm[p2].vp[I2_1][J2_1], V4 = model.csswm[p2].vp[I2_2][J2_2];

                alpha = model.alpha2D[I1][J1];
                beta = model.beta2D[I1][J1];

                uIP = model.interpolate(A1, A2, V1, V2, B);
                vIP = model.interpolate(A1, A2, V3, V4, B);

                alpha_B = model.interpolate(A1, A2, model.alpha2D[I2_1][J2_1], model.alpha2D[I2_2][J2_2], B);
                beta_B = model.interpolate(A1, A2, model.beta2D[I2_1][J2_1], model.beta2D[I2_2][J2_2], B);

                model.get_gLower(gLower, alpha, beta);
                model.get_IA(IA, p1, alpha, beta);
                model.get_A(A, p2, alpha_B, beta_B);
                model.get_gUpper(gUpper, alpha_B, beta_B);
                model.csswm[p1].up[I1][J1] = model.Cube2Cube_U_2(gLower, IA, A, gUpper, uIP, vIP);
                model.csswm[p1].vp[I1][J1] = model.Cube2Cube_V_2(gLower, IA, A, gUpper, uIP, vIP);
                
            }
            else {
                int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                int I2_1 = i2 == -1 ? reversed ? checkIP[NX-1-idx][0] : checkIP[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? checkIP[NY-1-idx][0] : checkIP[idx][0] : j2;
                int I2_2 = i2 == -1 ? reversed ? checkIP[NX-1-idx][1] : checkIP[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? checkIP[NY-1-idx][1] : checkIP[idx][1] : j2;

                B = model.csswm[p1].lon[I1][J1];
                A1 = model.csswm[p2].lon[I2_1][J2_1], A2 = model.csswm[p2].lon[I2_2][J2_2];
                V1 = model.csswm[p2].up[I2_1][J2_1], V2 = model.csswm[p2].up[I2_2][J2_2];
                V3 = model.csswm[p2].vp[I2_1][J2_1], V4 = model.csswm[p2].vp[I2_2][J2_2];

                if (A1 > A2 && (p1 == 0 || p2 == 0))  A2 += 2 * M_PI;
                if (A1 > B && B < A2) B += 2 * M_PI;
                
                alpha = model.alpha2D[I1][J1];
                beta = model.beta2D[I1][J1];

                uIP = model.interpolate(A1, A2, V1, V2, B);
                vIP = model.interpolate(A1, A2, V3, V4, B);

                alpha_B = model.interpolate(A1, A2, model.alpha2D[I2_1][J2_1], model.alpha2D[I2_2][J2_2], B);
                beta_B = model.interpolate(A1, A2, model.beta2D[I2_1][J2_1], model.beta2D[I2_2][J2_2], B);

                model.get_gLower(gLower, alpha, beta);
                model.get_IA(IA, p1, alpha, beta);
                model.get_A(A, p2, alpha_B, beta_B);
                model.get_gUpper(gUpper, alpha_B, beta_B);

                // std::cout << gUpper[0] << " " << gUpper[1] << " " << gUpper[2] << " " << gUpper[3] << " " << std::endl;
                model.csswm[p1].up[I1][J1] = model.Cube2Cube_U_2(gLower, IA, A, gUpper, uIP, vIP);
                model.csswm[p1].vp[I1][J1] = model.Cube2Cube_V_2(gLower, IA, A, gUpper, uIP, vIP);
            }
        }
    }
}
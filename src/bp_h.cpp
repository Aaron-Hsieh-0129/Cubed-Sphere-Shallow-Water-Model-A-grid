#include "constrcution.hpp"

void CSSWM::BP_h(CSSWM &model) {
    double B, A1, A2, V1, V2;
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
    for (int pp = 0; pp < 24; pp++) {
        p1 = match[pp][0], p2 = match[pp][1], i1 = match[pp][2], j1 = match[pp][3], i2 = match[pp][4], j2 = match[pp][5], reversed = match[pp][6], lonlat = match[pp][7];
        for (int idx = 0; idx < NX; idx++) {
            if (lonlat == 0) {
                int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                int I2_1 = i2 == -1 ? reversed ? checkIP[NX-1-idx][0] : checkIP[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? checkIP[NY-1-idx][0] : checkIP[idx][0] : j2;
                int I2_2 = i2 == -1 ? reversed ? checkIP[NX-1-idx][1] : checkIP[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? checkIP[NY-1-idx][1] : checkIP[idx][1] : j2;

                B = model.csswm[p1].lat[I1][J1];
                A1 = model.csswm[p2].lat[I2_1][J2_1], A2 = model.csswm[p2].lat[I2_2][J2_2];
                V1 = model.csswm[p2].hp[I2_1][J2_1], V2 = model.csswm[p2].hp[I2_2][J2_2];
                
                model.csswm[p1].hp[I1][J1] = interpolate(A1, A2, V1, V2, B);
            }
            else {
                int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                int I2_1 = i2 == -1 ? reversed ? checkIP[NX-1-idx][0] : checkIP[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? checkIP[NY-1-idx][0] : checkIP[idx][0] : j2;
                int I2_2 = i2 == -1 ? reversed ? checkIP[NX-1-idx][1] : checkIP[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? checkIP[NY-1-idx][1] : checkIP[idx][1] : j2;

                B = model.csswm[p1].lon[I1][J1];
                A1 = model.csswm[p2].lon[I2_1][J2_1], A2 = model.csswm[p2].lon[I2_2][J2_2];
                V1 = model.csswm[p2].hp[I2_1][J2_1], V2 = model.csswm[p2].hp[I2_2][J2_2];

                if (A1 > A2 && (p1 == 0 || p2 == 0))  A2 += 2 * M_PI;
                if (A1 > B && B < A2) B += 2 * M_PI;
                // std::cout << p1 << " " << p2 << " " << I1 << " " << J1 << " " << A1 << " " << A2 << std::endl;
                
                model.csswm[p1].hp[I1][J1] = interpolate(A1, A2, V1, V2, B);
            }
        }
    }
}

double CSSWM::interpolate(double A1, double A2, double V1, double V2, double B) {
    if (A1 == A2) {
        return V1;
    }

    if (!((A1 < B && B < A2) || (A1 > B && B > A2))) {
        std::cout << "Error at interpolation: " << A1 << " " << B << " " << A2 << std::endl;
    }
    
    return V1 + (V2-V1) * (B-A1) / (A2-A1);
}
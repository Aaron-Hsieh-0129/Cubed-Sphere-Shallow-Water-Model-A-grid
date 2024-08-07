#include "construction.hpp"

void CSSWM::BP_h(CSSWM &model) {
    #if defined(SecondOrderSpace)
    double B, A1, A2, V1, V2;

    int p1, p2, i1, j1, i2, j2, reversed, lonlat;
    for (int pp = 0; pp < 24; pp++) {
        p1 = model.match[pp][0], p2 = model.match[pp][1], i1 = model.match[pp][2], j1 = model.match[pp][3], i2 = model.match[pp][4], j2 = model.match[pp][5], reversed = model.match[pp][6], lonlat = model.match[pp][7];
        for (int idx = 0; idx < model.nx; idx++) {
            if (lonlat == 0) {
                int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                int I2_1 = i2 == -1 ? reversed ? model.checkIP[model.nx-1-idx][0] : model.checkIP[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP[model.ny-1-idx][0] : model.checkIP[idx][0] : j2;
                int I2_2 = i2 == -1 ? reversed ? model.checkIP[model.nx-1-idx][1] : model.checkIP[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP[model.ny-1-idx][1] : model.checkIP[idx][1] : j2;

                B = model.lat[p1][I1][J1];
                A1 = model.lat[p2][I2_1][J2_1], A2 = model.lat[p2][I2_2][J2_2];
                V1 = model.hp[p2][I2_1][J2_1], V2 = model.hp[p2][I2_2][J2_2];
                
                model.hp[p1][I1][J1] = interpolate(A1, A2, V1, V2, B);
            }
            else {
                int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                int I2_1 = i2 == -1 ? reversed ? model.checkIP[model.nx-1-idx][0] : model.checkIP[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP[model.ny-1-idx][0] : model.checkIP[idx][0] : j2;
                int I2_2 = i2 == -1 ? reversed ? model.checkIP[model.nx-1-idx][1] : model.checkIP[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP[model.ny-1-idx][1] : model.checkIP[idx][1] : j2;

                B = model.lon[p1][I1][J1];
                A1 = model.lon[p2][I2_1][J2_1], A2 = model.lon[p2][I2_2][J2_2];
                V1 = model.hp[p2][I2_1][J2_1], V2 = model.hp[p2][I2_2][J2_2];

                if (A1 > A2 && (p1 == 0 || p2 == 0))  A2 += 2 * M_PI;
                if (A1 > B && B < A2) B += 2 * M_PI;
                // std::cout << p1 << " " << p2 << " " << I1 << " " << J1 << " " << A1 << " " << A2 << std::endl;
                
                model.hp[p1][I1][J1] = interpolate(A1, A2, V1, V2, B);
            }
        }
    }

    #elif defined(FourthOrderSpace)
        double B, A1, A2, V1, V2;

        // Interpolate ouTTer
        int p1, p2, i1, j1, i2, j2, reversed, lonlat;
        for (int pp = 0; pp < 24; pp++) {
            p1 = model.match_ouTTer[pp][0], p2 = model.match_ouTTer[pp][1], i1 = model.match_ouTTer[pp][2], j1 = model.match_ouTTer[pp][3], i2 = model.match_ouTTer[pp][4], j2 = model.match_ouTTer[pp][5], reversed = model.match_ouTTer[pp][6], lonlat = model.match_ouTTer[pp][7];
            for (int idx = 0; idx < model.nx; idx++) {
                if (lonlat == 0) {
                    int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                    int I2_1 = i2 == -1 ? reversed ? model.checkIP_ouTTer[model.nx-1-idx][0] : model.checkIP_ouTTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP_ouTTer[model.ny-1-idx][0] : model.checkIP_ouTTer[idx][0] : j2;
                    int I2_2 = i2 == -1 ? reversed ? model.checkIP_ouTTer[model.nx-1-idx][1] : model.checkIP_ouTTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP_ouTTer[model.ny-1-idx][1] : model.checkIP_ouTTer[idx][1] : j2;

                    B = model.lat[p1][I1][J1];
                    A1 = model.lat[p2][I2_1][J2_1], A2 = model.lat[p2][I2_2][J2_2];
                    V1 = model.hp[p2][I2_1][J2_1], V2 = model.hp[p2][I2_2][J2_2];
                    
                    model.hp[p1][I1][J1] = interpolate(A1, A2, V1, V2, B);
                }
                else {
                    int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                    int I2_1 = i2 == -1 ? reversed ? model.checkIP_ouTTer[model.nx-1-idx][0] : model.checkIP_ouTTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP_ouTTer[model.ny-1-idx][0] : model.checkIP_ouTTer[idx][0] : j2;
                    int I2_2 = i2 == -1 ? reversed ? model.checkIP_ouTTer[model.nx-1-idx][1] : model.checkIP_ouTTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP_ouTTer[model.ny-1-idx][1] : model.checkIP_ouTTer[idx][1] : j2;

                    B = model.lon[p1][I1][J1];
                    A1 = model.lon[p2][I2_1][J2_1], A2 = model.lon[p2][I2_2][J2_2];
                    V1 = model.hp[p2][I2_1][J2_1], V2 = model.hp[p2][I2_2][J2_2];

                    if (A1 > A2 && (p1 == 0 || p2 == 0))  A2 += 2 * M_PI;
                    if (A1 > B && B < A2) B += 2 * M_PI;
                    // std::cout << p1 << " " << p2 << " " << I1 << " " << J1 << " " << A1 << " " << A2 << std::endl;
                    
                    model.hp[p1][I1][J1] = interpolate(A1, A2, V1, V2, B);
                }
            }
        }

        // Interpolate ouTer
        for (int pp = 0; pp < 24; pp++) {
            p1 = model.match_ouTer[pp][0], p2 = model.match_ouTer[pp][1], i1 = model.match_ouTer[pp][2], j1 = model.match_ouTer[pp][3], i2 = model.match_ouTer[pp][4], j2 = model.match_ouTer[pp][5], reversed = model.match_ouTer[pp][6], lonlat = model.match_ouTer[pp][7];
            for (int idx = 0; idx < model.nx; idx++) {
                if (lonlat == 0) {
                    int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                    int I2_1 = i2 == -1 ? reversed ? model.checkIP_ouTer[model.nx-1-idx][0] : model.checkIP_ouTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP_ouTer[model.ny-1-idx][0] : model.checkIP_ouTer[idx][0] : j2;
                    int I2_2 = i2 == -1 ? reversed ? model.checkIP_ouTer[model.nx-1-idx][1] : model.checkIP_ouTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP_ouTer[model.ny-1-idx][1] : model.checkIP_ouTer[idx][1] : j2;

                    B = model.lat[p1][I1][J1];
                    A1 = model.lat[p2][I2_1][J2_1], A2 = model.lat[p2][I2_2][J2_2];
                    V1 = model.hp[p2][I2_1][J2_1], V2 = model.hp[p2][I2_2][J2_2];
                    
                    model.hp[p1][I1][J1] = interpolate(A1, A2, V1, V2, B);
                }
                else {
                    int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                    int I2_1 = i2 == -1 ? reversed ? model.checkIP_ouTer[model.nx-1-idx][0] : model.checkIP_ouTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP_ouTer[model.ny-1-idx][0] : model.checkIP_ouTer[idx][0] : j2;
                    int I2_2 = i2 == -1 ? reversed ? model.checkIP_ouTer[model.nx-1-idx][1] : model.checkIP_ouTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP_ouTer[model.ny-1-idx][1] : model.checkIP_ouTer[idx][1] : j2;

                    B = model.lon[p1][I1][J1];
                    A1 = model.lon[p2][I2_1][J2_1], A2 = model.lon[p2][I2_2][J2_2];
                    V1 = model.hp[p2][I2_1][J2_1], V2 = model.hp[p2][I2_2][J2_2];

                    if (A1 > A2 && (p1 == 0 || p2 == 0))  A2 += 2 * M_PI;
                    if (A1 > B && B < A2) B += 2 * M_PI;
                    // std::cout << p1 << " " << p2 << " " << I1 << " " << J1 << " " << A1 << " " << A2 << std::endl;
                    
                    model.hp[p1][I1][J1] = interpolate(A1, A2, V1, V2, B);
                }
            }
        }
    #endif
}

#if defined(Mountain)
void CSSWM::BP_hs(CSSWM &model) {
    #if defined(SecondOrderSpace)
    double B, A1, A2, V1, V2;

    int p1, p2, i1, j1, i2, j2, reversed, lonlat;
    for (int pp = 0; pp < 24; pp++) {
        p1 = model.match[pp][0], p2 = model.match[pp][1], i1 = model.match[pp][2], j1 = model.match[pp][3], i2 = model.match[pp][4], j2 = model.match[pp][5], reversed = model.match[pp][6], lonlat = model.match[pp][7];
        for (int idx = 0; idx < model.nx; idx++) {
            if (lonlat == 0) {
                int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                int I2_1 = i2 == -1 ? reversed ? model.checkIP[model.nx-1-idx][0] : model.checkIP[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP[model.ny-1-idx][0] : model.checkIP[idx][0] : j2;
                int I2_2 = i2 == -1 ? reversed ? model.checkIP[model.nx-1-idx][1] : model.checkIP[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP[model.ny-1-idx][1] : model.checkIP[idx][1] : j2;

                B = model.lat[p1][I1][J1];
                A1 = model.lat[p2][I2_1][J2_1], A2 = model.lat[p2][I2_2][J2_2];
                V1 = model.hs[p2][I2_1][J2_1], V2 = model.hs[p2][I2_2][J2_2];
                
                model.hs[p1][I1][J1] = interpolate(A1, A2, V1, V2, B);
            }
            else {
                int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                int I2_1 = i2 == -1 ? reversed ? model.checkIP[model.nx-1-idx][0] : model.checkIP[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP[model.ny-1-idx][0] : model.checkIP[idx][0] : j2;
                int I2_2 = i2 == -1 ? reversed ? model.checkIP[model.nx-1-idx][1] : model.checkIP[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP[model.ny-1-idx][1] : model.checkIP[idx][1] : j2;

                B = model.lon[p1][I1][J1];
                A1 = model.lon[p2][I2_1][J2_1], A2 = model.lon[p2][I2_2][J2_2];
                V1 = model.hs[p2][I2_1][J2_1], V2 = model.hs[p2][I2_2][J2_2];

                if (A1 > A2 && (p1 == 0 || p2 == 0))  A2 += 2 * M_PI;
                if (A1 > B && B < A2) B += 2 * M_PI;
                // std::cout << p1 << " " << p2 << " " << I1 << " " << J1 << " " << A1 << " " << A2 << std::endl;
                
                model.hs[p1][I1][J1] = interpolate(A1, A2, V1, V2, B);
            }
        }
    }

    #elif defined(FourthOrderSpace)
        double B, A1, A2, V1, V2;

        // Interpolate ouTTer
        int p1, p2, i1, j1, i2, j2, reversed, lonlat;
        for (int pp = 0; pp < 24; pp++) {
            p1 = model.match_ouTTer[pp][0], p2 = model.match_ouTTer[pp][1], i1 = model.match_ouTTer[pp][2], j1 = model.match_ouTTer[pp][3], i2 = model.match_ouTTer[pp][4], j2 = model.match_ouTTer[pp][5], reversed = model.match_ouTTer[pp][6], lonlat = model.match_ouTTer[pp][7];
            for (int idx = 0; idx < model.nx; idx++) {
                if (lonlat == 0) {
                    int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                    int I2_1 = i2 == -1 ? reversed ? model.checkIP_ouTTer[model.nx-1-idx][0] : model.checkIP_ouTTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP_ouTTer[model.ny-1-idx][0] : model.checkIP_ouTTer[idx][0] : j2;
                    int I2_2 = i2 == -1 ? reversed ? model.checkIP_ouTTer[model.nx-1-idx][1] : model.checkIP_ouTTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP_ouTTer[model.ny-1-idx][1] : model.checkIP_ouTTer[idx][1] : j2;

                    B = model.lat[p1][I1][J1];
                    A1 = model.lat[p2][I2_1][J2_1], A2 = model.lat[p2][I2_2][J2_2];
                    V1 = model.hs[p2][I2_1][J2_1], V2 = model.hs[p2][I2_2][J2_2];
                    
                    model.hs[p1][I1][J1] = interpolate(A1, A2, V1, V2, B);
                }
                else {
                    int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                    int I2_1 = i2 == -1 ? reversed ? model.checkIP_ouTTer[model.nx-1-idx][0] : model.checkIP_ouTTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP_ouTTer[model.ny-1-idx][0] : model.checkIP_ouTTer[idx][0] : j2;
                    int I2_2 = i2 == -1 ? reversed ? model.checkIP_ouTTer[model.nx-1-idx][1] : model.checkIP_ouTTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP_ouTTer[model.ny-1-idx][1] : model.checkIP_ouTTer[idx][1] : j2;

                    B = model.lon[p1][I1][J1];
                    A1 = model.lon[p2][I2_1][J2_1], A2 = model.lon[p2][I2_2][J2_2];
                    V1 = model.hs[p2][I2_1][J2_1], V2 = model.hs[p2][I2_2][J2_2];

                    if (A1 > A2 && (p1 == 0 || p2 == 0))  A2 += 2 * M_PI;
                    if (A1 > B && B < A2) B += 2 * M_PI;
                    // std::cout << p1 << " " << p2 << " " << I1 << " " << J1 << " " << A1 << " " << A2 << std::endl;
                    
                    model.hs[p1][I1][J1] = interpolate(A1, A2, V1, V2, B);
                }
            }
        }

        // Interpolate ouTer
        for (int pp = 0; pp < 24; pp++) {
            p1 = model.match_ouTer[pp][0], p2 = model.match_ouTer[pp][1], i1 = model.match_ouTer[pp][2], j1 = model.match_ouTer[pp][3], i2 = model.match_ouTer[pp][4], j2 = model.match_ouTer[pp][5], reversed = model.match_ouTer[pp][6], lonlat = model.match_ouTer[pp][7];
            for (int idx = 0; idx < model.nx; idx++) {
                if (lonlat == 0) {
                    int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                    int I2_1 = i2 == -1 ? reversed ? model.checkIP_ouTer[model.nx-1-idx][0] : model.checkIP_ouTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP_ouTer[model.ny-1-idx][0] : model.checkIP_ouTer[idx][0] : j2;
                    int I2_2 = i2 == -1 ? reversed ? model.checkIP_ouTer[model.nx-1-idx][1] : model.checkIP_ouTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP_ouTer[model.ny-1-idx][1] : model.checkIP_ouTer[idx][1] : j2;

                    B = model.lat[p1][I1][J1];
                    A1 = model.lat[p2][I2_1][J2_1], A2 = model.lat[p2][I2_2][J2_2];
                    V1 = model.hs[p2][I2_1][J2_1], V2 = model.hs[p2][I2_2][J2_2];
                    
                    model.hs[p1][I1][J1] = interpolate(A1, A2, V1, V2, B);
                }
                else {
                    int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                    int I2_1 = i2 == -1 ? reversed ? model.checkIP_ouTer[model.nx-1-idx][0] : model.checkIP_ouTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP_ouTer[model.ny-1-idx][0] : model.checkIP_ouTer[idx][0] : j2;
                    int I2_2 = i2 == -1 ? reversed ? model.checkIP_ouTer[model.nx-1-idx][1] : model.checkIP_ouTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP_ouTer[model.ny-1-idx][1] : model.checkIP_ouTer[idx][1] : j2;

                    B = model.lon[p1][I1][J1];
                    A1 = model.lon[p2][I2_1][J2_1], A2 = model.lon[p2][I2_2][J2_2];
                    V1 = model.hs[p2][I2_1][J2_1], V2 = model.hs[p2][I2_2][J2_2];

                    if (A1 > A2 && (p1 == 0 || p2 == 0))  A2 += 2 * M_PI;
                    if (A1 > B && B < A2) B += 2 * M_PI;
                    // std::cout << p1 << " " << p2 << " " << I1 << " " << J1 << " " << A1 << " " << A2 << std::endl;
                    
                    model.hs[p1][I1][J1] = interpolate(A1, A2, V1, V2, B);
                }
            }
        }
    #endif
}
#endif

double CSSWM::interpolate(double A1, double A2, double V1, double V2, double B) {
    if (A1 == A2) {
        return V1;
    }

    if (!((A1 < B && B < A2) || (A1 > B && B > A2))) {
        std::cout << "Error at interpolation: " << A1 * 180. / M_PI << " " << B * 180. / M_PI << " " << A2 * 180. / M_PI << std::endl;
    }
    
    return V1 + (V2-V1) * (B-A1) / (A2-A1);
}
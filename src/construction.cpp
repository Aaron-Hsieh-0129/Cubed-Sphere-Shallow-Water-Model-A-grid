#include "construction.hpp"

CSSWM::patch::patch() {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            hp[i][j] = h[i][j] = hm[i][j] = FILLVALUE;
            up[i][j] = u[i][j] = um[i][j] = FILLVALUE;
            vp[i][j] = v[i][j] = vm[i][j] = FILLVALUE;

            lon[i][j] = lat[i][j] = FILLVALUE;

            x[i][j] = y[i][j] = FILLVALUE;
        }
    }
}

void CSSWM::Construct_gamma_sqrtG_GUpper(double alpha2D[NX][NY], double beta2D[NX][NY], double gamma[NX][NY], double sqrtG[NX][NY], double gUpper[NX][NY][4], double gLower[NX][NY][4]) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            gamma[i][j] = sqrt(1 + pow(tan(alpha2D[i][j]), 2) + pow(tan(beta2D[i][j]), 2));
            sqrtG[i][j] = 1. / (pow(gamma[i][j], 3) * pow(cos(alpha2D[i][j]), 2) * pow(cos(beta2D[i][j]), 2));

            gUpper[i][j][0] = pow((gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j])), 2) * (1 + pow(tan(beta2D[i][j]), 2));
            gUpper[i][j][1] = pow((gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j])), 2) * (tan(alpha2D[i][j]) * tan(beta2D[i][j]));
            gUpper[i][j][2] = gUpper[i][j][1];
            gUpper[i][j][3] = pow((gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j])), 2) * (1 + pow(tan(alpha2D[i][j]), 2));

            gLower[i][j][0] = 1. / (pow(gamma[i][j], 4) * pow(cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2)) * (1 + pow(tan(alpha2D[i][j]), 2));
            gLower[i][j][1] = 1. / (pow(gamma[i][j], 4) * pow(cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2)) * (-tan(alpha2D[i][j]) * tan(beta2D[i][j]));
            gLower[i][j][2] = gLower[i][j][1];
            gLower[i][j][3] = 1. / (pow(gamma[i][j], 4) * pow(cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2)) * (1 + pow(tan(beta2D[i][j]), 2));
        }
    }
    return;
}

void CSSWM::Construct_p0123_lonlat_xy_AIA(int p, double alpha2D[NX][NY], double beta2D[NX][NY], double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            // lon/lat
            lon[i][j] = alpha2D[i][j] + p * M_PI/2.;
            lat[i][j] = atan(tan(beta2D[i][j]) * cos(alpha2D[i][j]));

            csswm[p].lon_original[i][j] = lon[i][j];

            // x/y
            x[i][j] = RADIUS * (lon[i][j] - p * M_PI/2.);
            y[i][j] = RADIUS * atan(tan(lat[i][j]) / cos(lon[i][j] - p * M_PI/2.));

            // A/IA
            A[i][j][0] = 1. / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * gamma[i][j] * cos(beta2D[i][j]);
            A[i][j][1] = 0.;
            A[i][j][2] = 1. / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-tan(alpha2D[i][j]) * sin(beta2D[i][j]));
            A[i][j][3] = 1. / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) / cos(beta2D[i][j]);

            IA[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) / cos(beta2D[i][j]);
            IA[i][j][1] = 0.;
            IA[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (tan(alpha2D[i][j]) * sin(beta2D[i][j]));
            IA[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]));

            if (lon[i][j] < 0) lon[i][j] += 2 * M_PI;
        }
    }
}

void CSSWM::Construct_p4_lonlat_xy_AIA(int p, double alpha2D[NX][NY], double beta2D[NX][NY], double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            // lon/lat
            lon[i][j] = atan2(tan(alpha2D[i][j]), -tan(beta2D[i][j]));
            lat[i][j] = atan(1 / sqrt(pow(tan(alpha2D[i][j]), 2)+pow(tan(beta2D[i][j]), 2)));

            csswm[p].lon_original[i][j] = lon[i][j];

            // x/y
            x[i][j] = RADIUS * atan(sin(lon[i][j]) / tan(lat[i][j]));
            y[i][j] = RADIUS * atan(-cos(lon[i][j]) / tan(lat[i][j]));

            // A/AInverse
            A[i][j][0] = 1. / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(lon[i][j]));
            A[i][j][1] = 1. / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(lon[i][j]));
            A[i][j][2] = 1. / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(lon[i][j]));
            A[i][j][3] = 1. / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(lon[i][j]));

            IA[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(lon[i][j]));
            IA[i][j][1] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (-gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(lon[i][j]));
            IA[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(lon[i][j]));
            IA[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(lon[i][j]));

            if (lon[i][j] < 0) lon[i][j] += 2 * M_PI;
        }
    }
}

void CSSWM::Construct_p5_lonlat_xy_AIA(int p, double alpha2D[NX][NY], double beta2D[NX][NY], double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            // lon/lat
            lon[i][j] = atan2(tan(alpha2D[i][j]), tan(beta2D[i][j]));
            lat[i][j] = -atan(1 / sqrt(pow(tan(alpha2D[i][j]), 2)+pow(tan(beta2D[i][j]), 2)));

            csswm[p].lon_original[i][j] = lon[i][j];

            // x/y
            x[i][j] = RADIUS * atan(-sin(lon[i][j]) / tan(lat[i][j]));
            y[i][j] = RADIUS * atan(-cos(lon[i][j]) / tan(lat[i][j]));

            // A/AInverse
            A[i][j][0] = 1. / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(lon[i][j]));
            A[i][j][1] = 1. / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(lon[i][j]));
            A[i][j][2] = 1. / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(lon[i][j]));
            A[i][j][3] = 1. / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(lon[i][j]));

            IA[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(lon[i][j]));
            IA[i][j][1] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(lon[i][j]));
            IA[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (-cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(lon[i][j]));
            IA[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(lon[i][j]));

            if (lon[i][j] < 0) lon[i][j] += 2 * M_PI;
        }
    }
}

#if defined(SecondOrderSpace)
void initMatch_1point(int match[24][8]) {
    // Construct a array for dealing with interpolation between all patch
    // (p1, p2, i1, j1, i2, j2, reversed, LonLat: 1/0) 
    // Left, Right, Up, Down
    int tmp[24][8] = {
        {0, 3, 0, -1, NX-2, -1, 0, 0},  {0, 1, NX-1, -1, 1, -1, 0, 0},    {0, 4, -1, NY-1, -1, 1, 0, 1},     {0, 5, -1, 0, -1, NY-2, 0, 1},
        {1, 0, 0, -1, NX-2, -1, 0, 0},  {1, 2, NX-1, -1, 1, -1, 0, 0},    {1, 4, -1, NY-1, NX-2, -1, 0, 1},  {1, 5, -1, 0, NX-2, -1, 1, 1},
        {2, 1, 0, -1, NX-2, -1, 0, 0},  {2, 3, NX-1, -1, 1, -1, 0, 0},    {2, 4, -1, NY-1, -1, NY-2, 1, 1},  {2, 5, -1, 0, -1, 1, 1, 1},
        {3, 2, 0, -1, NX-2, -1, 0, 0},  {3, 0, NX-1, -1, 1, -1, 0, 0},    {3, 4, -1, NY-1, 1, -1, 1, 1},     {3, 5, -1, 0, 1, -1, 0, 1},
        {4, 3, 0, -1, -1, NY-2, 1, 1},  {4, 1, NX-1, -1, -1, NY-2, 0, 1}, {4, 2, -1, NY-1, -1, NY-2, 1, 1},  {4, 0, -1, 0, -1, NY-2, 0, 1},
        {5, 3, 0, -1, -1, 1, 0, 1},     {5, 1, NX-1, -1, -1, 1, 1, 1},    {5, 0, -1, NY-1, -1, 1, 0, 1},     {5, 2, -1, 0, -1, 1, 1, 1}
    };
    for (int i = 0; i < 24; i++) {
        for (int j = 0; j < 8; j++) {
            match[i][j] = tmp[i][j];
        }
    }
    return;
}
#elif defined(FourthOrderSpace)
void initMatch_2point(int match_ouTer[24][8], int match_ouTTer[24][8]) {
    // Construct a array for dealing with interpolation between all patch
    // (p1, p2, i1, j1, i2, j2, reversed, LonLat: 1/0) 
    // Left, Right, Up, Down
    int tmp_ouTTer[24][8] = {
        {0, 3, 0, -1, NX-4, -1, 0, 0},  {0, 1, NX-1, -1, 3, -1, 0, 0},    {0, 4, -1, NY-1, -1, 3, 0, 1},     {0, 5, -1, 0, -1, NY-4, 0, 1},
        {1, 0, 0, -1, NX-4, -1, 0, 0},  {1, 2, NX-1, -1, 3, -1, 0, 0},    {1, 4, -1, NY-1, NX-4, -1, 0, 1},  {1, 5, -1, 0, NX-4, -1, 1, 1},
        {2, 1, 0, -1, NX-4, -1, 0, 0},  {2, 3, NX-1, -1, 3, -1, 0, 0},    {2, 4, -1, NY-1, -1, NY-4, 1, 1},  {2, 5, -1, 0, -1, 3, 1, 1},
        {3, 2, 0, -1, NX-4, -1, 0, 0},  {3, 0, NX-1, -1, 3, -1, 0, 0},    {3, 4, -1, NY-1, 3, -1, 1, 1},     {3, 5, -1, 0, 3, -1, 0, 1},
        {4, 3, 0, -1, -1, NY-4, 1, 1},  {4, 1, NX-1, -1, -1, NY-4, 0, 1}, {4, 2, -1, NY-1, -1, NY-4, 1, 1},  {4, 0, -1, 0, -1, NY-4, 0, 1},
        {5, 3, 0, -1, -1, 3, 0, 1},     {5, 1, NX-1, -1, -1, 3, 1, 1},    {5, 0, -1, NY-1, -1, 3, 0, 1},     {5, 2, -1, 0, -1, 3, 1, 1}
    };

    int tmp_ouTer[24][8] = {
        {0, 3, 1, -1, NX-3, -1, 0, 0},  {0, 1, NX-2, -1, 2, -1, 0, 0},    {0, 4, -1, NY-2, -1, 2, 0, 1},     {0, 5, -1, 1, -1, NY-3, 0, 1},
        {1, 0, 1, -1, NX-3, -1, 0, 0},  {1, 2, NX-2, -1, 2, -1, 0, 0},    {1, 4, -1, NY-2, NX-3, -1, 0, 1},  {1, 5, -1, 1, NX-3, -1, 1, 1},
        {2, 1, 1, -1, NX-3, -1, 0, 0},  {2, 3, NX-2, -1, 2, -1, 0, 0},    {2, 4, -1, NY-2, -1, NY-3, 1, 1},  {2, 5, -1, 1, -1, 2, 1, 1},
        {3, 2, 1, -1, NX-3, -1, 0, 0},  {3, 0, NX-2, -1, 2, -1, 0, 0},    {3, 4, -1, NY-2, 2, -1, 1, 1},     {3, 5, -1, 1, 2, -1, 0, 1},
        {4, 3, 1, -1, -1, NY-3, 1, 1},  {4, 1, NX-2, -1, -1, NY-3, 0, 1}, {4, 2, -1, NY-2, -1, NY-3, 1, 1},  {4, 0, -1, 1, -1, NY-3, 0, 1},
        {5, 3, 1, -1, -1, 2, 0, 1},     {5, 1, NX-2, -1, -1, 2, 1, 1},    {5, 0, -1, NY-2, -1, 2, 0, 1},     {5, 2, -1, 1, -1, 2, 1, 1}
    };

    for (int i = 0; i < 24; i++) {
        for (int j = 0; j < 8; j++) {
            match_ouTTer[i][j] = tmp_ouTTer[i][j];
            match_ouTer[i][j] = tmp_ouTer[i][j];
        }
    }
    return;
}
#endif

CSSWM::CSSWM() {
    // Init new 1D array
    double *alpha = new double[NX], *beta = new double[NY];

    for (int i = 0; i < NX; i++) {
        #if defined(SecondOrderSpace)
            alpha[i] = -M_PI/4. + (M_PI/2.) / (NX-2) * (i-0.5);
        #elif defined(FourthOrderSpace)
            alpha[i] = -M_PI/4. + (M_PI/2.) / (NX-4) * (i-1.5);
        #endif
    }
    for (int j = 0; j < NY; j++) {
        #if defined(SecondOrderSpace)
            beta[j] = -M_PI/4. + (M_PI/2.) / (NX-2) * (j-0.5);
        #elif defined(FourthOrderSpace)
            beta[j] = -M_PI/4. + (M_PI/2.) / (NX-4) * (j-1.5);
        #endif
    }


    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            alpha2D[i][j] = alpha[i];
            beta2D[i][j] = beta[j];
        }
    }

    Construct_gamma_sqrtG_GUpper(alpha2D, beta2D, gamma, sqrtG, gUpper, gLower);

    for (int p = 0; p < 6; p++) {
        if (p == 4) {
            Construct_p4_lonlat_xy_AIA(p, alpha2D, beta2D, gamma, csswm[p].lon, csswm[p].lat, csswm[p].x, csswm[p].y, csswm[p].A, csswm[p].IA);
            continue;
        }
        if (p == 5) {
            Construct_p5_lonlat_xy_AIA(p, alpha2D, beta2D, gamma, csswm[p].lon, csswm[p].lat, csswm[p].x, csswm[p].y, csswm[p].A, csswm[p].IA);
            continue;
        }
        Construct_p0123_lonlat_xy_AIA(p, alpha2D, beta2D, gamma, csswm[p].lon, csswm[p].lat, csswm[p].x, csswm[p].y, csswm[p].A, csswm[p].IA);
    }

    #if defined(SecondOrderSpace)
        // Interpolation matrix construction
        int idx = 0, ipIdx = 0;
        double A1, A2, B;

        // Construct Interpolation 2D array filled with index (Note that all interpolation relation between every boundary is same)
        while (idx < NX && ipIdx < NX - 1) {
            B = csswm[0].lat[NX-1][idx];
            A1 = csswm[1].lat[1][ipIdx], A2 = csswm[1].lat[1][ipIdx+1];

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
    #elif defined(FourthOrderSpace)
        // For ouTTer
        // Interpolation matrix construction
        int idx = 0, ipIdx = 0;
        double A1, A2, B;
        // Construct Interpolation 2D array filled with index (Note that all interpolation relation between every boundary is same)
        while (idx < NX && ipIdx < NX - 1) {
            B = csswm[0].lat[NX-1][idx];
            A1 = csswm[1].lat[3][ipIdx], A2 = csswm[1].lat[3][ipIdx+1];

            if (A1 < B && B < A2) {
                checkIP_ouTTer[idx][0] = ipIdx;
                checkIP_ouTTer[idx][1] = ipIdx + 1;
                idx++;
            }
            else if (A1 == B) {
                checkIP_ouTTer[idx][0] = ipIdx;
                checkIP_ouTTer[idx][1] = ipIdx;
                idx++;
            }
            else if (A2 == B) {
                checkIP_ouTTer[idx][0] = ipIdx + 1;
                checkIP_ouTTer[idx][1] = ipIdx + 1;
                idx++;
            }
            else {
                ipIdx++;
            }
        } 

        // For ouTer
        // Interpolation matrix construction
        idx = 0, ipIdx = 0;
        A1 = 0, A2 = 0, B = 0;
        // Construct Interpolation 2D array filled with index (Note that all interpolation relation between every boundary is same)
        while (idx < NX && ipIdx < NX - 1) {
            B = csswm[0].lat[NX-2][idx];
            A1 = csswm[1].lat[2][ipIdx], A2 = csswm[1].lat[2][ipIdx+1];

            if (A1 < B && B < A2) {
                checkIP_ouTer[idx][0] = ipIdx;
                checkIP_ouTer[idx][1] = ipIdx + 1;
                idx++;
            }
            else if (A1 == B) {
                checkIP_ouTer[idx][0] = ipIdx;
                checkIP_ouTer[idx][1] = ipIdx;
                idx++;
            }
            else if (A2 == B) {
                checkIP_ouTer[idx][0] = ipIdx + 1;
                checkIP_ouTer[idx][1] = ipIdx + 1;
                idx++;
            }
            else {
                ipIdx++;
            }
        }
    #endif

    // Construct a array for dealing with interpolation between all patch
    // (p1, p2, i1, j1, i2, j2, reversed, LonLat: 1/0) 
    // Left, Right, Up, Down
    #if defined(SecondOrderSpace)
        initMatch_1point(match);
    #elif defined(FourthOrderSpace)
        initMatch_2point(match_ouTer, match_ouTTer);
    #endif

    // Construct patch to patch transformation matrix (Declare at transform.cpp)
    // Cube2Cube_matrix();

    delete[] alpha;
    delete[] beta;


    // define mountain
    #ifdef Mountain
        for (int p = 0; p < 6; p++) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    double r0 = M_PI / 9.;
                    double lonC = 3. * M_PI / 2., latC = M_PI / 6.;
                    double where = sqrt(pow(csswm[p].lon[i][j] - lonC, 2) + pow(csswm[p].lat[i][j] - latC, 2));
                    double r = r0 >  where ? where : r0;
                    double hs0 = 2000.;
                    csswm[p].hs[i][j] = hs0 * (1 - r / r0);
                }
            }
        }
    #endif
}

void CSSWM::get_gUpper(double ans[4], double alpha, double beta) {
    double gamma = sqrt(1 + pow(tan(alpha), 2) + pow(tan(beta), 2));

    ans[0] = pow((gamma * cos(alpha) * cos(beta)), 2) * (1 + pow(tan(beta), 2));
    ans[1] = pow((gamma * cos(alpha) * cos(beta)), 2) * (tan(alpha) * tan(beta));
    ans[2] = ans[1];
    ans[3] = pow((gamma* cos(alpha) * cos(beta)), 2) * (1 + pow(tan(alpha), 2));
}

void CSSWM::get_gLower(double ans[4], double alpha, double beta) {
    double gamma = sqrt(1 + pow(tan(alpha), 2) + pow(tan(beta), 2));
    
    ans[0] = 1. / (pow(gamma, 4) * pow(cos(alpha) * cos(beta), 2)) * (1 + pow(tan(alpha), 2));
    ans[1] = 1. / (pow(gamma, 4) * pow(cos(alpha) * cos(beta), 2)) * (-tan(alpha) * tan(beta));
    ans[2] = ans[1];
    ans[3] = 1. / (pow(gamma, 4) * pow(cos(alpha) * cos(beta), 2)) * (1 + pow(tan(beta), 2));
}

void CSSWM::get_A(double ans[4], int p, double alpha, double beta) {
    if (p == 0 || p == 1 || p == 2 || p == 3) {
        double gamma = sqrt(1 + pow(tan(alpha), 2) + pow(tan(beta), 2));

        ans[0] = 1. / (pow(gamma, 2) * cos(alpha) * cos(beta)) * gamma * cos(beta);
        ans[1] = 0.;
        ans[2] = 1. / (pow(gamma, 2) * cos(alpha) * cos(beta)) * (-tan(alpha) * sin(beta));
        ans[3] = 1. / (pow(gamma, 2) * cos(alpha) * cos(beta)) / cos(beta);
        return;
    }
    else if (p == 4) {
        double gamma = sqrt(1 + pow(tan(alpha), 2) + pow(tan(beta), 2));
        double lon = atan2(tan(alpha), -tan(beta));

        ans[0] = 1. / (pow(gamma, 2) * cos(alpha) * cos(beta)) * (gamma * cos(beta) / cos(alpha) * cos(lon));
        ans[1] = 1. / (pow(gamma, 2) * cos(alpha) * cos(beta)) * (gamma * cos(alpha) / cos(beta) * sin(lon));
        ans[2] = 1. / (pow(gamma, 2) * cos(alpha) * cos(beta)) * (-cos(beta) / cos(alpha) * sin(lon));
        ans[3] = 1. / (pow(gamma, 2) * cos(alpha) * cos(beta)) * (cos(alpha) / cos(beta) * cos(lon));
        return;
    }
    else {
        double gamma = sqrt(1 + pow(tan(alpha), 2) + pow(tan(beta), 2));
        double lon = atan2(tan(alpha), tan(beta));

        ans[0] = 1. / (pow(gamma, 2) * cos(alpha) * cos(beta)) * (gamma * cos(beta) / cos(alpha) * cos(lon));
        ans[1] = 1. / (pow(gamma, 2) * cos(alpha) * cos(beta)) * (-gamma * cos(alpha) / cos(beta) * sin(lon));
        ans[2] = 1. / (pow(gamma, 2) * cos(alpha) * cos(beta)) * (cos(beta) / cos(alpha) * sin(lon));
        ans[3] = 1. / (pow(gamma, 2) * cos(alpha) * cos(beta)) * (cos(alpha) / cos(beta) * cos(lon));
        return;
    }
}


void CSSWM::get_IA(double ans[4], int p, double alpha, double beta) {
    if (p == 0 || p == 1 || p == 2 || p == 3) {
        double gamma = sqrt(1 + pow(tan(alpha), 2) + pow(tan(beta), 2));

        ans[0] = gamma * cos(alpha) * cos(beta) / cos(beta);
        ans[1] = 0;
        ans[2] = gamma * cos(alpha) * cos(beta) * (tan(alpha) * sin(beta));
        ans[3] = gamma * cos(alpha) * cos(beta) * (gamma * cos(beta));
        return;
    }
    else if (p == 4) {
        double gamma = sqrt(1 + pow(tan(alpha), 2) + pow(tan(beta), 2));
        double lon = atan2(tan(alpha), -tan(beta));

        ans[0] = gamma * cos(alpha) * cos(beta) * (cos(alpha) / cos(beta) * cos(lon));
        ans[1] = gamma * cos(alpha) * cos(beta) * (-gamma * cos(alpha) / cos(beta) * sin(lon));
        ans[2] = gamma * cos(alpha) * cos(beta) * (cos(beta) / cos(alpha) * sin(lon));
        ans[3] = gamma * cos(alpha) * cos(beta) * (gamma * cos(beta) / cos(alpha) * cos(lon));
        return;
    }
    else {
        double gamma = sqrt(1 + pow(tan(alpha), 2) + pow(tan(beta), 2));
        double lon = atan2(tan(alpha), tan(beta));

        ans[0] = gamma * cos(alpha) * cos(beta) * (cos(alpha) / cos(beta) * cos(lon));
        ans[1] = gamma * cos(alpha) * cos(beta) * (gamma * cos(alpha) / cos(beta) * sin(lon));
        ans[2] = gamma * cos(alpha) * cos(beta) * (-cos(beta) / cos(alpha) * sin(lon));
        ans[3] = gamma * cos(alpha) * cos(beta) * (gamma * cos(beta) / cos(alpha) * cos(lon));
        return;
    }   
}
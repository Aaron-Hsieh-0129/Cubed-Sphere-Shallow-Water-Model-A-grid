#include "constrcution.hpp"

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

void CSSWM::Construct_gamma_sqrtG_GUpper(double **alpha2D, double **beta2D, double gamma[NX][NY], double sqrtG[NX][NY], double gUpper[NX][NY][4], double gLower[NX][NY][4]) {
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

void CSSWM::Construct_p0123_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            // lon/lat
            lon[i][j] = alpha2D[i][j] + p * M_PI/2.;
            lat[i][j] = atan(tan(beta2D[i][j]) * cos(alpha2D[i][j]));

            // x/y
            x[i][j] = RADIUS * (lon[i][j] - p * M_PI/2.);
            y[i][j] = RADIUS * atan(tan(lat[i][j]) / cos(lon[i][j] - p * M_PI/2.));

            // A/IA
            A[i][j][0] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * gamma[i][j] * cos(beta2D[i][j]);
            A[i][j][1] = 0;
            A[i][j][2] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-tan(alpha2D[i][j]) * sin(beta2D[i][j]));
            A[i][j][3] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) / cos(beta2D[i][j]);

            IA[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) / cos(beta2D[i][j]);
            IA[i][j][1] = 0;
            IA[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (tan(alpha2D[i][j]) * sin(beta2D[i][j]));
            IA[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]));

            if (lon[i][j] < 0) lon[i][j] += 2 * M_PI;
        }
    }
}

void CSSWM::Construct_p4_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            // lon/lat
            lon[i][j] = atan2(tan(alpha2D[i][j]), -tan(beta2D[i][j]));
            lat[i][j] = atan(1 / sqrt(pow(tan(alpha2D[i][j]), 2)+pow(tan(beta2D[i][j]), 2)));

            // x/y
            x[i][j] = RADIUS * atan(sin(lon[i][j]) / tan(lat[i][j]));
            y[i][j] = RADIUS * atan(-cos(lon[i][j]) / tan(lat[i][j]));

            // A/AInverse
            A[i][j][0] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(lon[i][j]));
            A[i][j][1] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(lon[i][j]));
            A[i][j][2] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(lon[i][j]));
            A[i][j][3] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(lon[i][j]));

            IA[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(lon[i][j]));
            IA[i][j][1] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (-gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(lon[i][j]));
            IA[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(lon[i][j]));
            IA[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(lon[i][j]));

            if (lon[i][j] < 0) lon[i][j] += 2 * M_PI;
        }
    }
}

void CSSWM::Construct_p5_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            // lon/lat
            lon[i][j] = atan2(tan(alpha2D[i][j]), tan(beta2D[i][j]));
            lat[i][j] = -atan(1 / sqrt(pow(tan(alpha2D[i][j]), 2)+pow(tan(beta2D[i][j]), 2)));

            // x/y
            x[i][j] = RADIUS * atan(-sin(lon[i][j]) / tan(lat[i][j]));
            y[i][j] = RADIUS * atan(-cos(lon[i][j]) / tan(lat[i][j]));

            // A/AInverse
            A[i][j][0] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(lon[i][j]));
            A[i][j][1] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(lon[i][j]));
            A[i][j][2] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(lon[i][j]));
            A[i][j][3] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(lon[i][j]));

            IA[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(lon[i][j]));
            IA[i][j][1] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(lon[i][j]));
            IA[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (-cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(lon[i][j]));
            IA[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(lon[i][j]));

            if (lon[i][j] < 0) lon[i][j] += 2 * M_PI;
        }
    }
}

CSSWM::CSSWM() {
    // Init new 1D array
    double *alpha = new double[NX], *beta = new double[NY];

    for (int i = 0; i < NX; i++) {
        alpha[i] = -M_PI/4. + (M_PI/2.) / (NX-4) * (i-1.5);
    }
    for (int j = 0; j < NY; j++) {
        beta[j] = -M_PI/4. + (M_PI/2.) / (NX-4) * (j-1.5);
    }

    // Init new 2D array
    double **alpha2D = new double *[NX], **beta2D = new double *[NY];
    for (int i = 0; i < NX; i++) {
        alpha2D[i] = new double[NY];
        beta2D[i] = new double [NY];
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

    // delete dynamic array
    delete[] alpha, delete[] beta;
    for (int i = 0; i < NX; i++) {
        delete[] alpha2D[i], delete[] beta2D[i];
    }
    delete[] alpha2D, delete[] beta2D;
}
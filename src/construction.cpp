#include "construction.hpp"

void CSSWM::Construct_gamma_sqrtG_GUpper(double **alpha2D, double **beta2D, double **gamma, double **sqrtG, double ***gUpper, double ***gLower) {
    for (int i = 0; i < CSSWM::nx; i++) {
        for (int j = 0; j < CSSWM::ny; j++) {
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

void CSSWM::Construct_p0123_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double **gamma, double **lon, double **lat, double **lon_original, double **x, double **y, double ***A, double ***IA) {
    for (int i = 0; i < CSSWM::nx; i++) {
        for (int j = 0; j < CSSWM::ny; j++) {
            // lon/lat
            lon[i][j] = alpha2D[i][j] + p * M_PI/2.;
            lat[i][j] = atan(tan(beta2D[i][j]) * cos(alpha2D[i][j]));

            lon_original[i][j] = lon[i][j];

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

void CSSWM::Construct_p4_lonlat_xy_AIA(double **alpha2D, double **beta2D, double **gamma, double **lon, double **lat, double **lon_original, double **x, double **y, double ***A, double ***IA) {
    for (int i = 0; i < CSSWM::nx; i++) {
        for (int j = 0; j < CSSWM::ny; j++) {
            // lon/lat
            lon[i][j] = atan2(tan(alpha2D[i][j]), -tan(beta2D[i][j]));
            lat[i][j] = atan(1 / sqrt(pow(tan(alpha2D[i][j]), 2)+pow(tan(beta2D[i][j]), 2)));

            lon_original[i][j] = lon[i][j];

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

void CSSWM::Construct_p5_lonlat_xy_AIA(double **alpha2D, double **beta2D, double **gamma, double **lon, double **lat, double **lon_original, double **x, double **y, double ***A, double ***IA) {
    for (int i = 0; i < CSSWM::nx; i++) {
        for (int j = 0; j < CSSWM::ny; j++) {
            // lon/lat
            lon[i][j] = atan2(tan(alpha2D[i][j]), tan(beta2D[i][j]));
            lat[i][j] = -atan(1 / sqrt(pow(tan(alpha2D[i][j]), 2)+pow(tan(beta2D[i][j]), 2)));

            lon_original[i][j] = lon[i][j];

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
void CSSWM::initMatch_1point(int **match) {
    // Construct a array for dealing with interpolation between all patch
    // (p1, p2, i1, j1, i2, j2, reversed, LonLat: 1/0) 
    // Left, Right, Up, Down
    int tmp[24][8] = {
        {0, 3, 0, -1, CSSWM::nx-2, -1, 0, 0},  {0, 1, CSSWM::nx-1, -1, 1, -1, 0, 0},    {0, 4, -1, CSSWM::ny-1, -1, 1, 0, 1},     {0, 5, -1, 0, -1, CSSWM::ny-2, 0, 1},
        {1, 0, 0, -1, CSSWM::nx-2, -1, 0, 0},  {1, 2, CSSWM::nx-1, -1, 1, -1, 0, 0},    {1, 4, -1, CSSWM::ny-1, CSSWM::nx-2, -1, 0, 1},  {1, 5, -1, 0, CSSWM::nx-2, -1, 1, 1},
        {2, 1, 0, -1, CSSWM::nx-2, -1, 0, 0},  {2, 3, CSSWM::nx-1, -1, 1, -1, 0, 0},    {2, 4, -1, CSSWM::ny-1, -1, CSSWM::ny-2, 1, 1},  {2, 5, -1, 0, -1, 1, 1, 1},
        {3, 2, 0, -1, CSSWM::nx-2, -1, 0, 0},  {3, 0, CSSWM::nx-1, -1, 1, -1, 0, 0},    {3, 4, -1, CSSWM::ny-1, 1, -1, 1, 1},     {3, 5, -1, 0, 1, -1, 0, 1},
        {4, 3, 0, -1, -1, CSSWM::ny-2, 1, 1},  {4, 1, CSSWM::nx-1, -1, -1, CSSWM::ny-2, 0, 1}, {4, 2, -1, CSSWM::ny-1, -1, CSSWM::ny-2, 1, 1},  {4, 0, -1, 0, -1, CSSWM::ny-2, 0, 1},
        {5, 3, 0, -1, -1, 1, 0, 1},     {5, 1, CSSWM::nx-1, -1, -1, 1, 1, 1},    {5, 0, -1, CSSWM::ny-1, -1, 1, 0, 1},     {5, 2, -1, 0, -1, 1, 1, 1}
    };
    for (int i = 0; i < 24; i++) {
        for (int j = 0; j < 8; j++) {
            match[i][j] = tmp[i][j];
        }
    }
    return;
}
#elif defined(FourthOrderSpace)
void CSSWM::initMatch_2point(int **match_ouTer, int **match_ouTTer) {
    // Construct a array for dealing with interpolation between all patch
    // (p1, p2, i1, j1, i2, j2, reversed, LonLat: 1/0) 
    // Left, Right, Up, Down
    int tmp_ouTTer[24][8] = {
        {0, 3, 0, -1, CSSWM::nx-4, -1, 0, 0},  {0, 1, CSSWM::nx-1, -1, 3, -1, 0, 0},    {0, 4, -1, CSSWM::ny-1, -1, 3, 0, 1},     {0, 5, -1, 0, -1, CSSWM::ny-4, 0, 1},
        {1, 0, 0, -1, CSSWM::nx-4, -1, 0, 0},  {1, 2, CSSWM::nx-1, -1, 3, -1, 0, 0},    {1, 4, -1, CSSWM::ny-1, CSSWM::nx-4, -1, 0, 1},  {1, 5, -1, 0, CSSWM::nx-4, -1, 1, 1},
        {2, 1, 0, -1, CSSWM::nx-4, -1, 0, 0},  {2, 3, CSSWM::nx-1, -1, 3, -1, 0, 0},    {2, 4, -1, CSSWM::ny-1, -1, CSSWM::ny-4, 1, 1},  {2, 5, -1, 0, -1, 3, 1, 1},
        {3, 2, 0, -1, CSSWM::nx-4, -1, 0, 0},  {3, 0, CSSWM::nx-1, -1, 3, -1, 0, 0},    {3, 4, -1, CSSWM::ny-1, 3, -1, 1, 1},     {3, 5, -1, 0, 3, -1, 0, 1},
        {4, 3, 0, -1, -1, CSSWM::ny-4, 1, 1},  {4, 1, CSSWM::nx-1, -1, -1, CSSWM::ny-4, 0, 1}, {4, 2, -1, CSSWM::ny-1, -1, CSSWM::ny-4, 1, 1},  {4, 0, -1, 0, -1, CSSWM::ny-4, 0, 1},
        {5, 3, 0, -1, -1, 3, 0, 1},     {5, 1, CSSWM::nx-1, -1, -1, 3, 1, 1},    {5, 0, -1, CSSWM::ny-1, -1, 3, 0, 1},     {5, 2, -1, 0, -1, 3, 1, 1}
    };

    int tmp_ouTer[24][8] = {
        {0, 3, 1, -1, CSSWM::nx-3, -1, 0, 0},  {0, 1, CSSWM::nx-2, -1, 2, -1, 0, 0},    {0, 4, -1, CSSWM::ny-2, -1, 2, 0, 1},     {0, 5, -1, 1, -1, CSSWM::ny-3, 0, 1},
        {1, 0, 1, -1, CSSWM::nx-3, -1, 0, 0},  {1, 2, CSSWM::nx-2, -1, 2, -1, 0, 0},    {1, 4, -1, CSSWM::ny-2, CSSWM::nx-3, -1, 0, 1},  {1, 5, -1, 1, CSSWM::nx-3, -1, 1, 1},
        {2, 1, 1, -1, CSSWM::nx-3, -1, 0, 0},  {2, 3, CSSWM::nx-2, -1, 2, -1, 0, 0},    {2, 4, -1, CSSWM::ny-2, -1, CSSWM::ny-3, 1, 1},  {2, 5, -1, 1, -1, 2, 1, 1},
        {3, 2, 1, -1, CSSWM::nx-3, -1, 0, 0},  {3, 0, CSSWM::nx-2, -1, 2, -1, 0, 0},    {3, 4, -1, CSSWM::ny-2, 2, -1, 1, 1},     {3, 5, -1, 1, 2, -1, 0, 1},
        {4, 3, 1, -1, -1, CSSWM::ny-3, 1, 1},  {4, 1, CSSWM::nx-2, -1, -1, CSSWM::ny-3, 0, 1}, {4, 2, -1, CSSWM::ny-2, -1, CSSWM::ny-3, 1, 1},  {4, 0, -1, 1, -1, CSSWM::ny-3, 0, 1},
        {5, 3, 1, -1, -1, 2, 0, 1},     {5, 1, CSSWM::nx-2, -1, -1, 2, 1, 1},    {5, 0, -1, CSSWM::ny-2, -1, 2, 0, 1},     {5, 2, -1, 1, -1, 2, 1, 1}
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
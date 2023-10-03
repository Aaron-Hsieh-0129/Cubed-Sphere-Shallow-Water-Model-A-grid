#include "initialField.hpp"

void Init::Init2d(CSSWM & model) {
    for (int p = 0; p < 6; p++) {
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                #ifdef ConvergenceRate
                    model.csswm[p].hp[i][j] = ConvergenceRateH(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]);
                    model.csswm[p].up[i][j] = (model.gLower[i][j][0] * model.csswm[p].IA[i][j][0] + model.gLower[i][j][1] * model.csswm[p].IA[i][j][2]) * SteadyGeostrophyU(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]) + 
                                              (model.gLower[i][j][0] * model.csswm[p].IA[i][j][1] + model.gLower[i][j][1] * model.csswm[p].IA[i][j][3]) * SteadyGeostrophyV(model.csswm[p].lon[i][j]);
                    model.csswm[p].vp[i][j] = (model.gLower[i][j][2] * model.csswm[p].IA[i][j][0] + model.gLower[i][j][3] * model.csswm[p].IA[i][j][2]) * SteadyGeostrophyU(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]) + 
                                              (model.gLower[i][j][2] * model.csswm[p].IA[i][j][1] + model.gLower[i][j][3] * model.csswm[p].IA[i][j][3]) * SteadyGeostrophyV(model.csswm[p].lon[i][j]);
                #endif

                #if defined(Advection)
                    model.csswm[p].hp[i][j] = AdvectionH(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]);
                    model.csswm[p].up[i][j] = (model.gLower[i][j][0] * model.csswm[p].IA[i][j][0] + model.gLower[i][j][1] * model.csswm[p].IA[i][j][2]) * AdvectionU(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]) + 
                                              (model.gLower[i][j][0] * model.csswm[p].IA[i][j][1] + model.gLower[i][j][1] * model.csswm[p].IA[i][j][3]) * AdvectionV(model.csswm[p].lon[i][j]);
                    model.csswm[p].vp[i][j] = (model.gLower[i][j][2] * model.csswm[p].IA[i][j][0] + model.gLower[i][j][3] * model.csswm[p].IA[i][j][2]) * AdvectionU(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]) + 
                                              (model.gLower[i][j][2] * model.csswm[p].IA[i][j][1] + model.gLower[i][j][3] * model.csswm[p].IA[i][j][3]) * AdvectionV(model.csswm[p].lon[i][j]);
                #endif

                #if defined(Jung)
                    model.csswm[p].hp[i][j] = JungH(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]);
                    double mult[2][2];
                    model.matrixMul(model.gLower[i][j], model.csswm[p].IA[i][j], mult);
                    model.csswm[p].up[i][j] = mult[0][0] * JungU(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]) + 
                                              mult[0][1] * JungV(model.csswm[p].lon[i][j]);

                    model.matrixMul(model.gLower[i][j], model.csswm[p].IA[i][j], mult);
                    model.csswm[p].vp[i][j] = mult[1][0] * JungU(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]) + 
                                              mult[1][1] * JungV(model.csswm[p].lon[i][j]);
                #endif

                #ifdef DeformationalFlow
                    model.csswm[p].hp[i][j] = DeformationalFlowH(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]);
                #endif

                #ifdef GravityWave
                    model.csswm[p].hp[i][j] = Gravity(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]);
                    model.csswm[p].up[i][j] = 0.;
                    model.csswm[p].vp[i][j] = 0.;
                #endif

                #ifdef SteadyGeostrophy
                    model.csswm[p].hp[i][j] = SteadyGeostrophyH(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]);
                    model.csswm[p].up[i][j] = (model.gLower[i][j][0] * model.csswm[p].IA[i][j][0] + model.gLower[i][j][1] * model.csswm[p].IA[i][j][2]) * SteadyGeostrophyU(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]) + 
                                              (model.gLower[i][j][0] * model.csswm[p].IA[i][j][1] + model.gLower[i][j][1] * model.csswm[p].IA[i][j][3]) * SteadyGeostrophyV(model.csswm[p].lon[i][j]);
                    model.csswm[p].vp[i][j] = (model.gLower[i][j][2] * model.csswm[p].IA[i][j][0] + model.gLower[i][j][3] * model.csswm[p].IA[i][j][2]) * SteadyGeostrophyU(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]) + 
                                              (model.gLower[i][j][2] * model.csswm[p].IA[i][j][1] + model.gLower[i][j][3] * model.csswm[p].IA[i][j][3]) * SteadyGeostrophyV(model.csswm[p].lon[i][j]);
                #endif

                #ifdef Barotropic
                    // model.csswm[p].hp[i][j] = BarotropicHPrime(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]);
                    model.csswm[p].hp[i][j] = BarotropicH(model.csswm[p].lat[i][j]) + BarotropicHPrime(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]);
                    model.csswm[p].up[i][j] = (model.gLower[i][j][0] * model.csswm[p].IA[i][j][0] + model.gLower[i][j][1] * model.csswm[p].IA[i][j][2]) * BarotropicU(model.csswm[p].lat[i][j]) + 
                                              (model.gLower[i][j][0] * model.csswm[p].IA[i][j][1] + model.gLower[i][j][1] * model.csswm[p].IA[i][j][3]) * 0;
                    model.csswm[p].vp[i][j] = (model.gLower[i][j][2] * model.csswm[p].IA[i][j][0] + model.gLower[i][j][3] * model.csswm[p].IA[i][j][2]) * BarotropicU(model.csswm[p].lat[i][j]) + 
                                              (model.gLower[i][j][2] * model.csswm[p].IA[i][j][1] + model.gLower[i][j][3] * model.csswm[p].IA[i][j][3]) * 0;
                #endif

                #ifdef Mountain
                    model.csswm[p].hp[i][j] = MountainH(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]);
                    model.csswm[p].up[i][j] = (model.gLower[i][j][0] * model.csswm[p].IA[i][j][0] + model.gLower[i][j][1] * model.csswm[p].IA[i][j][2]) * MountainU(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]) + 
                                              (model.gLower[i][j][0] * model.csswm[p].IA[i][j][1] + model.gLower[i][j][1] * model.csswm[p].IA[i][j][3]) * MountainV(model.csswm[p].lon[i][j]);
                    model.csswm[p].vp[i][j] = (model.gLower[i][j][2] * model.csswm[p].IA[i][j][0] + model.gLower[i][j][3] * model.csswm[p].IA[i][j][2]) * MountainU(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]) + 
                                              (model.gLower[i][j][2] * model.csswm[p].IA[i][j][1] + model.gLower[i][j][3] * model.csswm[p].IA[i][j][3]) * MountainV(model.csswm[p].lon[i][j]); 
                #endif

                #ifdef RossbyHaurwitz
                    model.csswm[p].hp[i][j] = RossbyHaurwitzH(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]);
                    model.csswm[p].up[i][j] = (model.gLower[i][j][0] * model.csswm[p].IA[i][j][0] + model.gLower[i][j][1] * model.csswm[p].IA[i][j][2]) * RossbyHaurwitzU(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]) + 
                                              (model.gLower[i][j][0] * model.csswm[p].IA[i][j][1] + model.gLower[i][j][1] * model.csswm[p].IA[i][j][3]) * RossbyHaurwitzV(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]);
                    model.csswm[p].vp[i][j] = (model.gLower[i][j][2] * model.csswm[p].IA[i][j][0] + model.gLower[i][j][3] * model.csswm[p].IA[i][j][2]) * RossbyHaurwitzU(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]) + 
                                              (model.gLower[i][j][2] * model.csswm[p].IA[i][j][1] + model.gLower[i][j][3] * model.csswm[p].IA[i][j][3]) * RossbyHaurwitzV(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]); 
                #endif

                #ifdef EquatorialWave
                    model.csswm[p].hp[i][j] = 10000. + EquatorialWave(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]);
                    model.csswm[p].up[i][j] = 0.;
                    model.csswm[p].vp[i][j] = 0.;
                #endif
            }
        }
    }

    #if defined(Mountain)
        model.BP_hs(model);
    #endif


    model.BP_h(model);
    #ifndef Advection
        // model.BP_wind_convert(model);
        // model.BP_wind_interpolation(model);
        model.BP_wind_interpolation2(model);
    #endif
    
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                model.csswm[p].hm[i][j] = model.csswm[p].hp[i][j]; 
                model.csswm[p].um[i][j] = model.csswm[p].up[i][j]; 
                model.csswm[p].vm[i][j] = model.csswm[p].vp[i][j];   

                model.csswm[p].h[i][j] = model.csswm[p].hp[i][j]; 
                model.csswm[p].u[i][j] = model.csswm[p].up[i][j]; 
                model.csswm[p].v[i][j] = model.csswm[p].vp[i][j];  
            }
        }
    }
}



double Init::ConvergenceRateH(double lon, double lat) {
    int p = 4, k = 4;
    double lonP = atan((cos(lat) * sin(lon)) / (cos(lat)*cos(lon)*cos(ALPHA0) + sin(lat)*sin(ALPHA0)));
    double latP = asin(sin(lat) * cos(ALPHA0) - cos(lat) * cos(lon) * sin(ALPHA0));
    return pow(cos(latP), p) * sin(k * lonP);
}

double Init::DeformationalFlowH(double lon, double lat) {
    int rho0 = 3, gamma = 5;
    double lonP = atan((cos(lat) * sin(lon)) / (cos(lat)*cos(lon)*cos(ALPHA0) + sin(lat)*sin(ALPHA0)));
    return 1 - tanh(rho0 / gamma * sin(lonP));

}

double Init::JungH(double lon, double lat) {
    double h0 = 1000;
    double lonC = 0., latC = 0.;
    double rd = RADIUS * acos(sin(latC) * sin(lat) + cos(latC) * cos(lat) * cos(lon-lonC));
    double r0 = RADIUS / 3.;
    if (rd < r0) return h0 / 2. * (1 + cos(M_PI * rd / r0));
    else return 0.;
}

double Init::JungU(double lon, double lat) {
    double u0 = 2 * M_PI * RADIUS / (12. * 86400);
    double u = u0 * (cos(ALPHA0) * cos(lat) + sin(ALPHA0) * sin(lon) * sin(lat));
    return u;
}

double Init::JungV(double lon) {
    double u0 = 2 * M_PI * RADIUS / (12. * 86400);
    double v = u0 * sin(ALPHA0) * cos(lon);
    return v;
}

double Init::AdvectionH(double lon, double lat) {
    double h0 = 1000.;
    double lonC = 3. * M_PI / 2., latC = 0.;
    double rd = RADIUS * acos(sin(latC) * sin(lat) + cos(latC) * cos(lat) * cos(lon-lonC));
    double r0 = RADIUS / 3.;
    if (rd < r0) return h0 / 2. * (1 + cos(M_PI * rd / r0));
    else return 0.;
}

double Init::AdvectionU(double lon, double lat) {
    double u0 = 2 * M_PI * RADIUS / (12. * 86400);
    double u = u0 * (cos(ALPHA0) * cos(lat) + sin(ALPHA0) * cos(lon) * sin(lat));
    return u;
}

double Init::AdvectionV(double lon) {
    double u0 = 2 * M_PI * RADIUS / (12. * 86400);
    double v = - u0 * sin(ALPHA0) * sin(lon);
    return v;
}

double Init::Gravity(double lon, double lat) {
    double h0 = 1000;
    double lonC = 0., latC = 0.;
    double rd = RADIUS * acos(sin(latC) * sin(lat) + cos(latC) * cos(lat) * cos(lon-lonC));
    double r0 = RADIUS / 3.;
    if (rd < r0) return h0 + h0 / 2. * (1 + cos(M_PI * rd / r0));
    else return h0;
}

double Init::SteadyGeostrophyH(double lon, double lat) {
    double h0 = 2.94E4 / GRAVITY;
    double u0 = 2 * M_PI * RADIUS / (12. * 86400);
    return h0 - (RADIUS * OMEGA * u0 + u0 * u0 / 2.) * pow(-cos(lon) * cos(lat) * sin(ALPHA0) + sin(lat) * cos(ALPHA0), 2) / GRAVITY;
}

double Init::SteadyGeostrophyU(double lon, double lat) {
    double u0 = 2 * M_PI * RADIUS / (12. * 86400);
    return u0 * (cos(ALPHA0) * cos(lat) + sin(ALPHA0) * cos(lon) * sin(lat));
}

double Init::SteadyGeostrophyV(double lon) {
    double u0 = 2 * M_PI * RADIUS / (12. * 86400);
    return -u0 * sin(ALPHA0) * sin(lon);
}

double Init::MountainH(double lon, double lat) {
    double h0 = 5960;
    double u0 = 20;
    return h0 - (RADIUS * OMEGA * u0 + u0 * u0 / 2.) * pow(-cos(lon) * cos(lat) * sin(ALPHA0) + sin(lat) * cos(ALPHA0), 2) / GRAVITY;
}

double Init::MountainU(double lon, double lat) {
    double u0 = 20;
    double u = u0 * (cos(ALPHA0) * cos(lat) + sin(ALPHA0) * cos(lon) * sin(lat));
    return u;
}

double Init::MountainV(double lon) {
    double u0 = 20;
    double v = - u0 * sin(ALPHA0) * sin(lon);
    return v;
}

double Init::BarotropicH(double lat) {
    double h0 = 10000.;
    return h0 - simpson(-M_PI / 2, lat) / GRAVITY;
}


double Init::BarotropicHPrime(double lon, double lat) {
    double alpha = 1. / 3., beta = 1. / 15., theta2 = M_PI / 4., hHat = 120;
    if (-M_PI < lon && lon < M_PI) return hHat * cos(lat) * exp(-pow(lon / alpha, 2) - pow((theta2 - lat) / beta, 2));
    else return 0;
}

double Init::BarotropicU(double lat) {
    double theta0 = M_PI / 7.;
    double theta1 = M_PI / 2. - theta0;
    double en = exp(-4. / pow(theta1 - theta0, 2));
    double umax = 80.;
    
    if (lat <= theta0) return 0.;
    else if (theta0 < lat && lat < theta1) return umax / en * exp(1. / ((lat - theta0) * (lat - theta1)));
    else return 0.;
}

double Init::func(double x) {
    double f = 2 * OMEGA * sin(x);
    return RADIUS * BarotropicU(x) * (f + tan(x) / RADIUS * BarotropicU(x));
}

double Init::simpson(double a, double b) {
	double I2n = 0, h = b - a;
	double T2n = h * (func(a) + func(b)) / 2.;
	double In = T2n;
	double Tn;
	for (int n = 1; abs(I2n-In) > 1e-5; n+=n, h/=2.0) {
		In = I2n;
		Tn = T2n;
		double sigma = 0.;
		for (int k = 0; k < n; k++) {
			sigma += func(a+(k+0.5) * h);
		}
		T2n = (Tn + h * sigma) / 2.;
		I2n = (4 * T2n - Tn) / 3.;
	}
	return I2n;
}

double Init::RossbyHaurwitzH(double lon, double lat) {
    double h0 = 8E3;
    double omega = 7.848E-6, K = 7.848E-6;
    double R = 4.;
    double c = cos(lat);

    double A = omega / 2. * (2. * OMEGA + omega) * c*c + 
               0.25 * K*K * pow(c, 2*R) * ((R+1)*c*2 + (2*R*R - R - 2) - 2*R*R*pow(c, -2));
    double B = (2 * (OMEGA + omega) * K) / ((R+1) * (R+2)) * pow(c, R) * ((R*R+2*R+2) - pow((R+1)*c, 2));
    double C = 1. / 4. * K*K * pow(c, 2*R) * ((R+1)*c*c - R+2);

    return (h0 + (RADIUS*RADIUS*A + RADIUS*RADIUS*B*cos(R*lon) + RADIUS*RADIUS*C*cos(2*R*lon)) / GRAVITY);
}

double Init::RossbyHaurwitzU(double lon, double lat) {
    double R = 4.;
    double omega = 7.848E-6, K = 7.848E-6;
    double c = cos(lat);
    return (RADIUS * omega * cos(lat) + RADIUS * K * pow(c, R-1) * (R*pow(sin(lat), 2) - pow(cos(lat), 2)) * cos(RADIUS*lon));
}

double Init::RossbyHaurwitzV(double lon, double lat) {
    double R = 4.;
    double K = 7.848E-6;
    double c = cos(lat);
    return (-RADIUS * K * R * pow(c, R-1) * sin(lat) * sin(RADIUS*lon));
}

double Init::EquatorialWaveH(double lon, double lat) {
    double h0 = 100;
    double lonC = 0., latC = 0.;
    double rd = RADIUS * acos(sin(latC) * sin(lat) + cos(latC) * cos(lat) * cos(lon-lonC));
    double r0 = RADIUS / 4.;
    if (rd < r0) return -h0 / 2. * (1 + cos(M_PI * rd / r0));
    else return h0;
}
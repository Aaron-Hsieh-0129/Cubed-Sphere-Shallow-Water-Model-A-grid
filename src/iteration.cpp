#include "iteration.hpp"

void Iteration::ph_pt(CSSWM &model) {
    for (int p = 0; p < 6; p++) {
        double psqrtGHU_px = 0, psqrtGHU_py = 0, dx_for_h = 0, dy_for_h = 0;
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                dx_for_h = model.csswm[p].x[i+1][j] - model.csswm[p].x[i-1][j];
                dy_for_h = model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j-1];

                psqrtGHU_px = (1. / (model.sqrtG[i][j] * dx_for_h)) * 
                              ((model.sqrtG[i+1][j] * model.csswm[p].h[i+1][j] * (model.gUpper[i+1][j][0] * model.csswm[p].u[i+1][j] + model.gUpper[i+1][j][1] * model.csswm[p].v[i+1][j])) - 
                               (model.sqrtG[i-1][j] * model.csswm[p].h[i-1][j] * (model.gUpper[i-1][j][0] * model.csswm[p].u[i-1][j] + model.gUpper[i-1][j][1] * model.csswm[p].v[i-1][j])));

                psqrtGHU_py = (1. / (model.sqrtG[i][j] * dy_for_h)) * 
                              ((model.sqrtG[i][j+1] * model.csswm[p].h[i][j+1] * (model.gUpper[i][j+1][2] * model.csswm[p].u[i][j+1] + model.gUpper[i][j+1][3] * model.csswm[p].v[i][j+1])) - 
                               (model.sqrtG[i][j-1] * model.csswm[p].h[i][j-1] * (model.gUpper[i][j-1][2] * model.csswm[p].u[i][j-1] + model.gUpper[i][j-1][3] * model.csswm[p].v[i][j-1])));
            
                model.csswm[p].hp[i][j] = model.csswm[p].hm[i][j] + D2T * (-psqrtGHU_px - psqrtGHU_py);
                
                #ifdef DIFFUSION
                    model.csswm[p].hp[i][j] += D2T * KX * (model.csswm[p].hm[i+1][j] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i-1][j]) / pow(dx_for_h, 2) + 
                                               D2T * KY * (model.csswm[p].hm[i][j+1] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i][j-1]) / pow(dy_for_h, 2);
                #endif
            }
        }
    }
    return;
}

void Iteration::pu_pt(CSSWM &model) {
    double dx_for_u = 0, dy_for_u = 0;
    double pgH_px = 0, pU2_px = 0, pUV_px = 0, pV2_px = 0, rotationU = 0;
    double f;
    #ifdef Mountain
        double pgHs_px = 0.;
    #endif

    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                dx_for_u = model.csswm[p].x[i+1][j] - model.csswm[p].x[i-1][j];
                dy_for_u = model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j-1];

                #ifdef SteadyGeostrophy
                    f = 2 * OMEGA * (-cos(model.csswm[p].lon_original[i][j] * cos(model.csswm[p].lat[i][j] * sin(ALPHA0) + sin(model.csswm[p].lat[i][j] * cos(ALPHA0)))));
                #elif defined(Barotropic) || defined(RossbyHaurwitz)
                    f = 2 * OMEGA * sin(model.csswm[p].lat[i][j]);
                #else
                    f = 0;
                #endif
                
                pgH_px = GRAVITY / dx_for_u * (model.csswm[p].h[i+1][j] - model.csswm[p].h[i-1][j]);

                pU2_px = 0.5 / dx_for_u * 
                        (model.gUpper[i+1][j][0] * pow(model.csswm[p].u[i+1][j], 2) - 
                         model.gUpper[i-1][j][0] * pow(model.csswm[p].u[i-1][j], 2));

                pUV_px = 1 / dx_for_u * 
                        (model.gUpper[i+1][j][1] * model.csswm[p].u[i+1][j] * model.csswm[p].v[i+1][j] - 
                         model.gUpper[i-1][j][1] * model.csswm[p].u[i-1][j] * model.csswm[p].v[i-1][j]);

                pV2_px = 0.5 / dx_for_u * 
                        (model.gUpper[i+1][j][3] * pow(model.csswm[p].v[i+1][j], 2) -
                         model.gUpper[i-1][j][3] * pow(model.csswm[p].v[i-1][j], 2));

                rotationU = (((model.csswm[p].v[i+1][j] - model.csswm[p].v[i-1][j]) / dx_for_u) - 
                             ((model.csswm[p].u[i][j+1] - model.csswm[p].u[i][j-1]) / dy_for_u) + model.sqrtG[i][j] * f) * 
                            (model.gUpper[i][j][2] * model.csswm[p].u[i][j] + model.gUpper[i][j][3] * model.csswm[p].v[i][j]);
            

                #ifdef Mountain
                    pgHs_px = GRAVITY / dx_for_u * (model.csswm[p].hs[i+1][j] - model.csswm[p].hs[i-1][j]);
                    model.csswm[p].up[i][j] = model.csswm[p].um[i][j] + D2T * (-pgH_px - pgHs_px - pU2_px - pUV_px - pV2_px + rotationU);
                #else
                    model.csswm[p].up[i][j] = model.csswm[p].um[i][j] + D2T * (-pgH_px - pU2_px - pUV_px - pV2_px + rotationU);
                #endif

                #ifdef DIFFUSION
                    model.csswm[p].up[i][j] += D2T * KX * (model.csswm[p].um[i+1][j] - 2. * model.csswm[p].um[i][j] + model.csswm[p].um[i-1][j]) / pow(dx_for_u, 2) + 
                                               D2T * KY * (model.csswm[p].um[i][j+1] - 2. * model.csswm[p].um[i][j] + model.csswm[p].um[i][j-1]) / pow(dy_for_u, 2);
                #endif
            }
        }
    }
    return;
}

void Iteration::pv_pt(CSSWM &model) {
    double dx_for_v = 0, dy_for_v = 0;
    double pgH_py = 0, pU2_py = 0, pUV_py = 0, pV2_py = 0, rotationV = 0;
    double f;
    #ifdef Mountain
        double pgHs_py = 0.;
    #endif

    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                dx_for_v = model.csswm[p].x[i+1][j] - model.csswm[p].x[i-1][j];
                dy_for_v = model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j-1];

                #if defined(SteadyGeostrophy)
                    f = 2 * OMEGA * (-cos(model.csswm[p].lon_original[i][j] * cos(model.csswm[p].lat[i][j] * sin(ALPHA0) + sin(model.csswm[p].lat[i][j] * cos(ALPHA0)))));
                #elif defined(Barotropic) || defined(RossbyHaurwitz)
                    f = 2 * OMEGA * sin(model.csswm[p].lat[i][j]);
                #else
                    f = 0;
                #endif

                pgH_py = GRAVITY / dy_for_v * (model.csswm[p].h[i][j+1] - model.csswm[p].h[i][j-1]);

                pU2_py = 0.5 / dy_for_v * 
                        (model.gUpper[i][j+1][0] * pow(model.csswm[p].u[i][j+1], 2) - 
                         model.gUpper[i][j-1][0] * pow(model.csswm[p].u[i][j-1], 2));

                pUV_py = 1 / dy_for_v * 
                        (model.gUpper[i][j+1][1] * model.csswm[p].u[i][j+1] * model.csswm[p].v[i][j+1] - 
                         model.gUpper[i][j-1][1] * model.csswm[p].u[i][j-1] * model.csswm[p].v[i][j-1]);

                pV2_py = 0.5 / dy_for_v * 
                        (model.gUpper[i][j+1][3] * pow(model.csswm[p].v[i][j+1], 2) -
                         model.gUpper[i][j-1][3] * pow(model.csswm[p].v[i][j-1], 2));

                rotationV = (((model.csswm[p].v[i+1][j] - model.csswm[p].v[i-1][j]) / dx_for_v) - 
                             ((model.csswm[p].u[i][j+1] - model.csswm[p].u[i][j-1]) / dy_for_v) + model.sqrtG[i][j] * f) * 
                            (model.gUpper[i][j][0] * model.csswm[p].u[i][j] + model.gUpper[i][j][1] * model.csswm[p].v[i][j]);

                
                #ifdef Mountain
                    pgHs_py = GRAVITY / dy_for_v * (model.csswm[p].hs[i][j+1] - model.csswm[p].hs[i][j-1]);;
                    model.csswm[p].vp[i][j] = model.csswm[p].vm[i][j] + D2T * (-pgH_py - pgHs_py - pU2_py - pUV_py - pV2_py - rotationV);
                #else
                    model.csswm[p].vp[i][j] = model.csswm[p].vm[i][j] + D2T * (-pgH_py - pU2_py - pUV_py - pV2_py - rotationV);
                #endif

                #ifdef DIFFUSION
                    model.csswm[p].vp[i][j] += D2T * KX * (model.csswm[p].vm[i+1][j] - 2. * model.csswm[p].vm[i][j] + model.csswm[p].vm[i-1][j]) / pow(dx_for_v, 2) + 
                                               D2T * KY * (model.csswm[p].vm[i][j+1] - 2. * model.csswm[p].vm[i][j] + model.csswm[p].vm[i][j-1]) / pow(dy_for_v, 2);
                #endif
            }
        }
    }
    return;
}

void Iteration::leap_frog(CSSWM &model) {
    Outputs::output_parameter(model);
    int n = 0;
    double timenow = 0.;
    double temp = TIMEEND / DT;
    int nmax = (int) temp;

    while (n < nmax) {
        std::cout << n << std::endl;

        if (n % OUTPUTINTERVAL == 0) {
            Outputs::output_h(n, model);
            Outputs::output_u(n, model);
            Outputs::output_v(n, model);
        }

        n++;
        timenow = n * DT;

        // Integrate
        ph_pt(model);
        #ifndef ADVECTION
            pu_pt(model);
            pv_pt(model);
        #endif

        // Boundary exchange and interpolation
        model.BP_h(model);
        #ifndef ADVECTION
            // model.BP_wind_convert(model);
            // model.BP_wind_interpolation(model);
            model.BP_wind_interpolation2(model);
        #endif

        // Time filter
        #ifdef TIMEFILTER
            for (int p = 0; p < 6; p++) {
                for (int i = 0; i < NX; i++) {
                    for (int j = 0; j < NY; j++) {
                        model.csswm[p].h[i][j] += TIMETS * (model.csswm[p].hp[i][j] - 2 * model.csswm[p].h[i][j] + model.csswm[p].hm[i][j]);
                        model.csswm[p].u[i][j] += TIMETS * (model.csswm[p].up[i][j] - 2 * model.csswm[p].u[i][j] + model.csswm[p].um[i][j]);
                        model.csswm[p].v[i][j] += TIMETS * (model.csswm[p].vp[i][j] - 2 * model.csswm[p].v[i][j] + model.csswm[p].vm[i][j]);
                    }
                }
            }
        #endif

        // next step
        for (int p = 0; p < 6; p++) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    model.csswm[p].hm[i][j] = model.csswm[p].h[i][j];
                    model.csswm[p].h[i][j] = model.csswm[p].hp[i][j];

                    model.csswm[p].um[i][j] = model.csswm[p].u[i][j];
                    model.csswm[p].u[i][j] = model.csswm[p].up[i][j];

                    model.csswm[p].vm[i][j] = model.csswm[p].v[i][j];
                    model.csswm[p].v[i][j] = model.csswm[p].vp[i][j];
                }
            }
        }
    }
    
    return;
}
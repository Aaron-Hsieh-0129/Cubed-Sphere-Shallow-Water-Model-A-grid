#include "construction.hpp"

#if defined(SecondOrderSpace)
void CSSWM::Iteration::ph_pt_2(CSSWM &model) {
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

                #if defined(EquatorialWave)
                    if (model.status_add_forcing == true) model.csswm[p].hp[i][j] = model.csswm[p].hm[i][j] + D2T * (-psqrtGHU_px - psqrtGHU_py + model.csswm[p].h_forcing[i][j]);
                    else model.csswm[p].hp[i][j] = model.csswm[p].hm[i][j] + D2T * (-psqrtGHU_px - psqrtGHU_py);
                #else
                    model.csswm[p].hp[i][j] = model.csswm[p].hm[i][j] + D2T * (-psqrtGHU_px - psqrtGHU_py);
                #endif
                
                #ifdef DIFFUSION
                    model.csswm[p].hp[i][j] += D2T * KX * (model.csswm[p].hm[i+1][j] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i-1][j]) / pow(dx_for_h, 2) + 
                                               D2T * KY * (model.csswm[p].hm[i][j+1] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i][j-1]) / pow(dy_for_h, 2);
                #endif
            }
        }
    }
    return;
}

void CSSWM::Iteration::pu_pt_2(CSSWM &model) {
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

                #if defined(SteadyGeostrophy) || defined(Mountain)
                    f = 2 * OMEGA * (-cos(model.csswm[p].lon[i][j]) * cos(model.csswm[p].lat[i][j]) * sin(ALPHA0) + sin(model.csswm[p].lat[i][j]) * cos(ALPHA0));
                #elif defined(Barotropic) || defined(RossbyHaurwitz)
                    f = 2 * OMEGA * sin(model.csswm[p].lat[i][j]);
                #elif defined(EquatorialWave)
                    double f0 = 0;
                    double beta = 2.5 * 10E-5;
                    f = f0 + beta * abs(model.csswm[p].lat[i][j] * 180 / M_PI) * (111000);
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

void CSSWM::Iteration::pv_pt_2(CSSWM &model) {
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

                #if defined(SteadyGeostrophy) || defined(Mountain)
                    f = 2 * OMEGA * (-cos(model.csswm[p].lon[i][j]) * cos(model.csswm[p].lat[i][j]) * sin(ALPHA0) + sin(model.csswm[p].lat[i][j]) * cos(ALPHA0));
                #elif defined(Barotropic) || defined(RossbyHaurwitz)
                    f = 2 * OMEGA * sin(model.csswm[p].lat[i][j]);
                #elif defined(EquatorialWave)
                    double f0 = 0;
                    double beta = 2.5 * 10E-5;
                    f = f0 + beta * abs(model.csswm[p].lat[i][j] * 180 / M_PI) * (111000);
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
#elif defined(FourthOrderSpace)
void CSSWM::Iteration::ph_pt_4(CSSWM &model) {
    for (int p = 0; p < 6; p++) {
        double psqrtGHU_px = 0, psqrtGHU_py = 0, dx_for_h = 0, dy_for_h = 0;
        for (int i = 2; i < NX-2; i++) {
            for (int j = 2; j < NY-2; j++) {
                dx_for_h = 0.5 * (model.csswm[p].x[i+1][j] - model.csswm[p].x[i-1][j]);
                dy_for_h = 0.5 * (model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j-1]);

                psqrtGHU_px = (1. / (model.sqrtG[i][j] * 12. * dx_for_h)) * 
                              (-1.*(model.sqrtG[i+2][j] * model.csswm[p].h[i+2][j] * (model.gUpper[i+2][j][0] * model.csswm[p].u[i+2][j] + model.gUpper[i+2][j][1] * model.csswm[p].v[i+2][j]))
                               +8.*(model.sqrtG[i+1][j] * model.csswm[p].h[i+1][j] * (model.gUpper[i+1][j][0] * model.csswm[p].u[i+1][j] + model.gUpper[i+1][j][1] * model.csswm[p].v[i+1][j]))
                               -8.*(model.sqrtG[i-1][j] * model.csswm[p].h[i-1][j] * (model.gUpper[i-1][j][0] * model.csswm[p].u[i-1][j] + model.gUpper[i-1][j][1] * model.csswm[p].v[i-1][j]))
                               +1.*(model.sqrtG[i-2][j] * model.csswm[p].h[i-2][j] * (model.gUpper[i-2][j][0] * model.csswm[p].u[i-2][j] + model.gUpper[i-2][j][1] * model.csswm[p].v[i-2][j])));

                psqrtGHU_py = (1. / (model.sqrtG[i][j] * 12. * dy_for_h)) * 
                              (-1.*(model.sqrtG[i][j+2] * model.csswm[p].h[i][j+2] * (model.gUpper[i][j+2][2] * model.csswm[p].u[i][j+2] + model.gUpper[i][j+2][3] * model.csswm[p].v[i][j+2])) 
                               +8.*(model.sqrtG[i][j+1] * model.csswm[p].h[i][j+1] * (model.gUpper[i][j+1][2] * model.csswm[p].u[i][j+1] + model.gUpper[i][j+1][3] * model.csswm[p].v[i][j+1]))
                               -8.*(model.sqrtG[i][j-1] * model.csswm[p].h[i][j-1] * (model.gUpper[i][j-1][2] * model.csswm[p].u[i][j-1] + model.gUpper[i][j-1][3] * model.csswm[p].v[i][j-1]))
                               +1.*(model.sqrtG[i][j-2] * model.csswm[p].h[i][j-2] * (model.gUpper[i][j-2][2] * model.csswm[p].u[i][j-2] + model.gUpper[i][j-2][3] * model.csswm[p].v[i][j-2])));
            
                #if defined(AB2Time)
                    model.csswm[p].dh[i][j][(model.step+1)%2] = (-psqrtGHU_px - psqrtGHU_py);
                    #if defined(EquatorialWave)
                        if (model.status_add_forcing == true) model.csswm[p].dh[i][j][(model.step+1)%2] += model.h_forcing[p][i][j];
                    #endif
                    if (model.step == 0) model.csswm[p].dh[i][j][0] = model.csswm[p].dh[i][j][1];
                    model.csswm[p].hp[i][j] = model.csswm[p].h[i][j] + 1.5*DT*model.csswm[p].dh[i][j][(model.step+1)%2] - 0.5*DT*model.csswm[p].dh[i][j][model.step%2];
                #else
                    #if defined(EquatorialWave)
                        if (model.status_add_forcing == true) model.csswm[p].hp[i][j] = model.csswm[p].hm[i][j] + D2T * (-psqrtGHU_px - psqrtGHU_py + model.csswm[p].h_forcing[i][j]);
                        else model.csswm[p].hp[i][j] = model.csswm[p].hm[i][j] + D2T * (-psqrtGHU_px - psqrtGHU_py);
                    #else
                        model.csswm[p].hp[i][j] = model.csswm[p].hm[i][j] + D2T * (-psqrtGHU_px - psqrtGHU_py);
                    #endif
                #endif

                // if ((p == 1 || p == 0) && j == 2) {
                //     model.csswm[p].hp[i][j] += 10. * -std::sin(i * M_PI / 2.);
                // }
            }
        }
    }
    return;
}

void CSSWM::Iteration::pu_pt_4(CSSWM &model) {
    double dx_for_u = 0, dy_for_u = 0;
    double pgH_px = 0, pU2_px = 0, pUV_px = 0, pV2_px = 0, rotationU = 0;
    double f;
    #ifdef Mountain
        double pgHs_px = 0.;
    #endif

    for (int p = 0; p < 6; p++) {
        for (int i = 2; i < NX-2; i++) {
            for (int j = 2; j < NY-2; j++) {
                dx_for_u = 0.5 * (model.csswm[p].x[i+1][j] - model.csswm[p].x[i-1][j]);
                dy_for_u = 0.5 * (model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j-1]);

                #if defined(SteadyGeostrophy) || defined(Mountain)
                    f = 2 * OMEGA * (-cos(model.csswm[p].lon[i][j]) * cos(model.csswm[p].lat[i][j]) * sin(ALPHA0) + sin(model.csswm[p].lat[i][j]) * cos(ALPHA0));
                #elif defined(Barotropic) || defined(RossbyHaurwitz)
                    f = 2 * OMEGA * sin(model.csswm[p].lat[i][j]);
                #elif defined(EquatorialWave)
                    double f0 = 0;
                    double beta = 2.5 * 10E-11;
                    f = f0 + beta * model.csswm[p].lat[i][j] * 180. / M_PI * (111000.);
                #else
                    f = 0;
                #endif
                
                pgH_px = GRAVITY / (12.*dx_for_u) * (-1*model.csswm[p].h[i+2][j] + 8*model.csswm[p].h[i+1][j] - 8*model.csswm[p].h[i-1][j] + 1*model.csswm[p].h[i-2][j]);

                pU2_px = 0.5 / (12.*dx_for_u) * 
                        (-1.*(model.gUpper[i+2][j][0] * pow(model.csswm[p].u[i+2][j], 2))
                         +8.*(model.gUpper[i+1][j][0] * pow(model.csswm[p].u[i+1][j], 2))
                         -8.*(model.gUpper[i-1][j][0] * pow(model.csswm[p].u[i-1][j], 2))
                         +1.*(model.gUpper[i-2][j][0] * pow(model.csswm[p].u[i-2][j], 2)));

                pUV_px = 1./ (12.*dx_for_u) * 
                        (-1.*(model.gUpper[i+2][j][1] * model.csswm[p].u[i+2][j] * model.csswm[p].v[i+2][j])
                         +8.*(model.gUpper[i+1][j][1] * model.csswm[p].u[i+1][j] * model.csswm[p].v[i+1][j])
                         -8.*(model.gUpper[i-1][j][1] * model.csswm[p].u[i-1][j] * model.csswm[p].v[i-1][j])
                         +1.*(model.gUpper[i-2][j][1] * model.csswm[p].u[i-2][j] * model.csswm[p].v[i-2][j]));

                pV2_px = 0.5 / (12.*dx_for_u) * 
                        (-1.*(model.gUpper[i+2][j][3] * pow(model.csswm[p].v[i+2][j], 2))
                         +8.*(model.gUpper[i+1][j][3] * pow(model.csswm[p].v[i+1][j], 2))
                         -8.*(model.gUpper[i-1][j][3] * pow(model.csswm[p].v[i-1][j], 2))
                         +1.*(model.gUpper[i-2][j][3] * pow(model.csswm[p].v[i-2][j], 2)));

                rotationU = (((-1.*model.csswm[p].v[i+2][j] + 8.*model.csswm[p].v[i+1][j] - 8.*model.csswm[p].v[i-1][j] + 1.*model.csswm[p].v[i-2][j]) / (12.*dx_for_u)) - 
                             ((-1.*model.csswm[p].u[i][j+2] + 8.*model.csswm[p].u[i][j+1] - 8.*model.csswm[p].u[i][j-1] + 1.*model.csswm[p].u[i][j-2]) / (12.*dy_for_u)) 
                             + model.sqrtG[i][j] * f) 
                             * (model.gUpper[i][j][2] * model.csswm[p].u[i][j] + model.gUpper[i][j][3] * model.csswm[p].v[i][j]);
            

                #if defined(AB2Time)
                    model.csswm[p].du[i][j][(model.step+1)%2] = (-pgH_px - pU2_px - pUV_px - pV2_px + rotationU);
                    #if defined(Mountain)
                        pgHs_px = model.gravity / (12.*dx_for_u) * (-1.*model.hs[p][i+2][j] + 8.*model.hs[p][i+1][j] - 8.*model.hs[p][i-1][j] + 1.*model.hs[p][i-2][j]);
                        model.csswm[p].du[i][j][(model.step+1)%2] += -pgHs_px;
                    #endif
                    if (model.step == 0) model.csswm[p].du[i][j][0] = model.csswm[p].du[i][j][1];
                    model.csswm[p].up[i][j] = model.csswm[p].u[i][j] + 1.5*DT*model.csswm[p].du[i][j][(model.step+1)%2] - 0.5*DT*model.csswm[p].du[i][j][model.step%2];
                #else
                    #if defined(Mountain)
                        pgHs_px = model.gravity / (12.*dx_for_u) * (-1.*model.hs[p][i+2][j] + 8.*model.hs[p][i+1][j] - 8.*model.hs[p][i-1][j] + 1.*model.hs[p][i-2][j]);
                        model.csswm[p].up[i][j] = model.um[p][i][j] + D2T * (-pgH_px - pgHs_px - pU2_px - pUV_px - pV2_px + rotationU);
                    #else
                        model.csswm[p].up[i][j] = model.csswm[p].um[i][j] + D2T * (-pgH_px - pU2_px - pUV_px - pV2_px + rotationU);
                    #endif
                #endif
            }
        }
    }
    return;
}

void CSSWM::Iteration::pv_pt_4(CSSWM &model) {
    double dx_for_v = 0, dy_for_v = 0;
    double pgH_py = 0, pU2_py = 0, pUV_py = 0, pV2_py = 0, rotationV = 0;
    double f;
    #ifdef Mountain
        double pgHs_py = 0.;
    #endif

    for (int p = 0; p < 6; p++) {
        for (int i = 2; i < NX-2; i++) {
            for (int j = 2; j < NY-2; j++) {
                dx_for_v = 0.5 * (model.csswm[p].x[i+1][j] - model.csswm[p].x[i-1][j]);
                dy_for_v = 0.5 * (model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j-1]);

                #if defined(SteadyGeostrophy) || defined(Mountain)
                    f = 2 * OMEGA * (-cos(model.csswm[p].lon[i][j]) * cos(model.csswm[p].lat[i][j]) * sin(ALPHA0) + sin(model.csswm[p].lat[i][j]) * cos(ALPHA0));
                #elif defined(Barotropic) || defined(RossbyHaurwitz)
                    f = 2 * OMEGA * sin(model.csswm[p].lat[i][j]);
                #elif defined(EquatorialWave)
                    double f0 = 0;
                    double beta = 2.5 * 10E-11;
                    f = f0 + beta * model.csswm[p].lat[i][j] * 180. / M_PI * (111000.);
                #else
                    f = 0;
                #endif

                pgH_py = GRAVITY / (12.*dy_for_v) * (-1.*model.csswm[p].h[i][j+2] + 8.*model.csswm[p].h[i][j+1] - 8.*model.csswm[p].h[i][j-1] + 1.*model.csswm[p].h[i][j-2]);

                pU2_py = 0.5 /(12.*dy_for_v) * 
                        (-1.*(model.gUpper[i][j+2][0] * pow(model.csswm[p].u[i][j+2], 2))
                         +8.*(model.gUpper[i][j+1][0] * pow(model.csswm[p].u[i][j+1], 2))
                         -8.*(model.gUpper[i][j-1][0] * pow(model.csswm[p].u[i][j-1], 2))
                         +1.*(model.gUpper[i][j-2][0] * pow(model.csswm[p].u[i][j-2], 2)));

                pUV_py = 1. / (12.*dy_for_v) * 
                        (-1.*(model.gUpper[i][j+2][1] * model.csswm[p].u[i][j+2] * model.csswm[p].v[i][j+2])
                         +8.*(model.gUpper[i][j+1][1] * model.csswm[p].u[i][j+1] * model.csswm[p].v[i][j+1])
                         -8.*(model.gUpper[i][j-1][1] * model.csswm[p].u[i][j-1] * model.csswm[p].v[i][j-1])
                         +1.*(model.gUpper[i][j-2][1] * model.csswm[p].u[i][j-2] * model.csswm[p].v[i][j-2]));

                pV2_py = 0.5 / (12.*dy_for_v) * 
                        (-1.*(model.gUpper[i][j+2][3] * pow(model.csswm[p].v[i][j+2], 2))
                         +8.*(model.gUpper[i][j+1][3] * pow(model.csswm[p].v[i][j+1], 2))
                         -8.*(model.gUpper[i][j-1][3] * pow(model.csswm[p].v[i][j-1], 2))
                         +1.*(model.gUpper[i][j-2][3] * pow(model.csswm[p].v[i][j-2], 2)));

                rotationV = (((-1.*model.csswm[p].v[i+2][j] + 8.*model.csswm[p].v[i+1][j] - 8.*model.csswm[p].v[i-1][j] + 1.*model.csswm[p].v[i-2][j]) / (12.*dx_for_v)) - 
                             ((-1.*model.csswm[p].u[i][j+2] + 8.*model.csswm[p].u[i][j+1] - 8.*model.csswm[p].u[i][j-1] + 1.*model.csswm[p].u[i][j-2]) / (12.*dy_for_v)) 
                             + model.sqrtG[i][j] * f) * 
                            (model.gUpper[i][j][0] * model.csswm[p].u[i][j] + model.gUpper[i][j][1] * model.csswm[p].v[i][j]);

                

                #if defined(AB2Time)
                    model.csswm[p].dv[i][j][(model.step+1)%2] = (-pgH_py - pU2_py - pUV_py - pV2_py - rotationV);
                    #if defined(Mountain)
                        pgHs_py = model.gravity / (12.*dy_for_u) * (-1.*model.hs[p][i+2][j] + 8.*model.hs[p][i+1][j] - 8.*model.hs[p][i-1][j] + 1.*model.hs[p][i-2][j]);
                        model.csswm[p].dv[i][j][(model.step+1)%2] += -pgHs_py;
                    #endif
                    if (model.step == 0) model.csswm[p].dv[i][j][0] = model.csswm[p].dv[i][j][1];
                    model.csswm[p].vp[i][j] = model.csswm[p].v[i][j] + 1.5*DT*model.csswm[p].dv[i][j][(model.step+1)%2] - 0.5*DT*model.csswm[p].dv[i][j][model.step%2];
                #else
                    #if defined(Mountain)
                        pgHs_py = model.gravity / (12.*dy_for_v) * (-1.*model.hs[p][i][j+2] + 8.*model.hs[p][i][j+1] - 8.*model.hs[p][i][j-1] + 1.*model.hs[p][i][j-2]);
                        model.csswm[p].vp[i][j] = model.vm[p][i][j] + D2T * (-pgH_py - pgHs_py - pU2_py - pUV_py - pV2_py - rotationV);
                    #else
                        model.csswm[p].vp[i][j] = model.csswm[p].vm[i][j] + D2T * (-pgH_py - pU2_py - pUV_py - pV2_py - rotationV);
                    #endif
                #endif
            }
        }
    }
    return;
}

#endif

void CSSWM::Iteration::nextTimeStep(CSSWM &model) {
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
    return;
}

void CSSWM::Iteration::TimeMarching(CSSWM &model) {
    Outputs::create_all_directory();
    #ifdef NCOUTPUT
        Outputs::grid_nc(model);
    #endif
    #ifdef TXTOUTPUT
        Outputs::grid(model);
    #endif
    model.step = 0;
    // double timenow = 0.;
    double temp = TIMEEND / DT;
    int nmax = (int) temp;

    while (model.step < nmax) {
        std::cout << model.step << std::endl;

        if (model.step % OUTPUTINTERVAL == 0) {
            #ifdef TXTOUTPUT
                Outputs::h(model.step, model);
                Outputs::u(model.step, model);
                Outputs::v(model.step, model);
            #endif

            #ifdef NCOUTPUT
                Outputs::huv_nc(model.step, model);
            #endif
        }

        model.step++;
        #if defined(EquatorialWave)
            if (model.step * DT >= ADDFORCINGTIME) model.status_add_forcing = false;
            else model.status_add_forcing = true;
        #endif


        // Integrate
        #if defined(SecondOrderSpace)
            ph_pt_2(model);
            #ifndef Advection
                pu_pt_2(model);
                pv_pt_2(model);
            #endif

            // Boundary exchange and interpolation
            model.BP_h(model);
            #ifndef Advection
                // model.BP_wind_convert(model);
                // model.BP_wind_interpolation(model);
                model.BP_wind_interpolation2(model);
            #endif
        #elif defined(FourthOrderSpace)
            ph_pt_4(model);
            #ifndef Advection
                pu_pt_4(model);
                pv_pt_4(model);
            #endif

            // Boundary exchange and interpolation
            model.BP_h(model);
            #ifndef Advection
                // model.BP_wind_convert(model);
                // model.BP_wind_interpolation(model);
                model.BP_wind_interpolation2(model);
            #endif
        #endif

        #if defined(DIFFUSION)
            CSSWM::NumericalProcess::DiffusionAll(model);
        #endif

        // // Time filter
        #if defined(TIMEFILTER) && !defined(AB2Time)
            CSSWM::NumericalProcess::timeFilterAll(model);
        #endif

        // Add forcing
        // Coupling time: 600
        // if (model.step == 3) {
        //     model.csswm[1].hp[46][47] += 6.518550756778495;
        //     model.csswm[1].hp[47][47] += 6.518550756778495;
        //     model.csswm[1].hp[48][47] += 6.518550756778495;
        // }
        // if (model.step == 6) {
        //     model.csswm[1].hp[46][47] += 13.576829193507365;
        //     model.csswm[1].hp[47][47] += 13.594794503515004;
        //     model.csswm[1].hp[48][47] += 13.576805050761322;
        // }
        // if (model.step == 9) {
        //     model.csswm[1].hp[46][47] += 5.630701927977498;
        //     model.csswm[1].hp[47][47] += 5.655678397986776;
        //     model.csswm[1].hp[48][47] += 5.630617705381155;
        // }
        // if (model.step == 12) {
        //     model.csswm[1].hp[46][47] += 4.4994444012518215;
        //     model.csswm[1].hp[47][47] += 4.468312874076219;
        //     model.csswm[1].hp[48][47] += 4.499312003414161;
        // }
        // if (model.step == 15) {
        //     model.csswm[1].hp[46][47] += -1.179906883449803;
        //     model.csswm[1].hp[47][47] += -1.2728558726466872;
        //     model.csswm[1].hp[48][47] += -1.18004420300349;
        // }
        // if (model.step == 18) {
        //     model.csswm[1].hp[46][47] += 1.0899403250605246;
        //     model.csswm[1].hp[47][47] += 0.8626470096041885;
        //     model.csswm[1].hp[48][47] += 1.0902238598555414;
        // }
        // if (model.step == 21) {
        //     model.csswm[1].hp[46][47] += 1.7982856258804532;
        //     model.csswm[1].hp[47][47] += 1.3580435658677743;
        //     model.csswm[1].hp[48][47] += 1.7986760395087913;
        // }

        // Coupling time: 1200s
        // if (model.step == 6) {
        //     model.csswm[1].hp[46][47] += 20.31932687189692;
        //     model.csswm[1].hp[47][47] += 20.31932687189692;
        //     model.csswm[1].hp[48][47] += 20.31932687189692;
        // }
        // if (model.step == 12) {
        //     model.csswm[1].hp[46][47] += 10.944539267993605;
        //     model.csswm[1].hp[47][47] += 10.911059770220163;
        //     model.csswm[1].hp[48][47] += 10.944441891570023;
        // }
        // if (model.step == 18) {
        //     model.csswm[1].hp[46][47] += 1.028993711599469;
        //     model.csswm[1].hp[47][47] += 0.8432120151355775;
        //     model.csswm[1].hp[48][47] += 1.0291502518593916;
        // }
        // if (model.step == 24) {
        //     model.csswm[1].hp[46][47] += 3.272221580031328;
        //     model.csswm[1].hp[47][47] += 2.342710919521778;
        //     model.csswm[1].hp[48][47] += 3.27377919268838;
        // }

        // Coupling time: 1800s
        // if (model.step == 9) {
        //     model.csswm[1].hp[46][47] += 26.55747411457378;
        //     model.csswm[1].hp[47][47] += 26.55747411457378;
        //     model.csswm[1].hp[48][47] += 26.55747411457378;
        // }
        // if (model.step == 18) {
        //     model.csswm[1].hp[46][47] += 7.4260124174943485;
        //     model.csswm[1].hp[47][47] += 7.32227240750035;
        //     model.csswm[1].hp[48][47] += 7.426011678808209;
        // }
        // if (model.step == 27) {
        //     model.csswm[1].hp[46][47] += 1.342739204685131;
        //     model.csswm[1].hp[47][47] += 1.2110326880829234;
        //     model.csswm[1].hp[48][47] += 1.3429452419659356;
        // }

        // Coupling time: 2400s
        // if (model.step == 12) {
        //     model.csswm[1].hp[46][47] += 31.86735402;
        //     model.csswm[1].hp[47][47] += 31.86735402;
        //     model.csswm[1].hp[48][47] += 31.86735402;
        // }
        // if (model.step == 24) {
        //     model.csswm[1].hp[46][47] += 3.65462477;
        //     model.csswm[1].hp[47][47] += 3.51916444;
        //     model.csswm[1].hp[48][47] += 3.65460265;
        // }

        // Coupling time: 3000s
        // if (model.step == 15) {
        //     model.csswm[1].hp[46][47] += 31.813260145421737;
        //     model.csswm[1].hp[47][47] += 31.813260145421737;
        //     model.csswm[1].hp[48][47] += 31.813260145421737;
        // }
        // if (model.step == 30) {
        //     model.csswm[1].hp[46][47] += 2.7420342661043833;
        //     model.csswm[1].hp[47][47] += 2.6085480010151514;
        //     model.csswm[1].hp[48][47] += 2.742061058903346;
        // }

        // Coupling time: 3600s
        // if (model.step == 18) {
        //     model.csswm[1].hp[46][47] += 33.489597494404734;
        //     model.csswm[1].hp[47][47] += 33.489597494404734;
        //     model.csswm[1].hp[48][47] += 33.489597494404734;
        // }


        // next step
        CSSWM::Iteration::nextTimeStep(model);
    }
    
    return;
}

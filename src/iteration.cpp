#include "construction.hpp"

#if defined(SecondOrderSpace)
void CSSWM::Iteration::ph_pt_2(CSSWM &model) {
    for (int p = 0; p < 6; p++) {
        double psqrtGHU_px = 0, psqrtGHU_py = 0, dx_for_h = 0, dy_for_h = 0;
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                dx_for_h = model.x[p][i+1][j] - model.x[p][i-1][j];
                dy_for_h = model.y[p][i][j+1] - model.y[p][i][j-1];

                psqrtGHU_px = (1. / (model.sqrtG[i][j] * dx_for_h)) * 
                              ((model.sqrtG[i+1][j] * model.h[p][i+1][j] * (model.gUpper[i+1][j][0] * model..u[p][i+1][j] + model.gUpper[i+1][j][1] * model.v[p][i+1][j])) - 
                               (model.sqrtG[i-1][j] * model.h[p][i-1][j] * (model.gUpper[i-1][j][0] * model.u[p][i-1][j] + model.gUpper[i-1][j][1] * model.v[p][i-1][j])));

                psqrtGHU_py = (1. / (model.sqrtG[i][j] * dy_for_h)) * 
                              ((model.sqrtG[i][j+1] * model.h[p][i][j+1] * (model.gUpper[i][j+1][2] * model.u[p][i][j+1] + model.gUpper[i][j+1][3] * model.v[p][i][j+1])) - 
                               (model.sqrtG[i][j-1] * model.h[p][i][j-1] * (model.gUpper[i][j-1][2] * model.u[p][i][j-1] + model.gUpper[i][j-1][3] * model.v[p][i][j-1])));

                #if defined(EquatorialWave)
                    if (model.status_add_forcing == true) model.hp[p][i][j] = model.hm[p][i][j] + model.d2t * (-psqrtGHU_px - psqrtGHU_py + model.h_forcing[p][i][j]);
                    else model.hp[p][i][j] = model.hm[p][i][j] + model.d2t * (-psqrtGHU_px - psqrtGHU_py);
                #else
                    model.hp[p][i][j] = model.hm[p][i][j] + model.d2t * (-psqrtGHU_px - psqrtGHU_py);
                #endif
                
                #ifdef DIFFUSION
                    model.hp[p][i][j] += model.d2t * KX * (model.hm[p][i+1][j] - 2. * model.hm[p][i][j] + model.hm[p][i-1][j]) / pow(dx_for_h, 2) + 
                                               model.d2t * KY * (model.hm[p][i][j+1] - 2. * model.hm[p][i][j] + model.hm[p][i][j-1]) / pow(dy_for_h, 2);
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
                dx_for_u = model.x[p][i+1][j] - model.x[p][i-1][j];
                dy_for_u = model.y[p][i][j+1] - model.y[p][i][j-1];

                #if defined(SteadyGeostrophy) || defined(Mountain)
                    f = 2 * OMEGA * (-cos(model.lon[p][i][j]) * cos(model.lat[p][i][j]) * sin(ALPHA0) + sin(model.lat[p][i][j]) * cos(ALPHA0));
                #elif defined(Barotropic) || defined(RossbyHaurwitz)
                    f = 2 * OMEGA * sin(model.lat[p][i][j]);
                #elif defined(EquatorialWave)
                    double f0 = 0;
                    double beta = 2.5 * 10E-5;
                    f = f0 + beta * abs(model.lat[p][i][j] * 180 / M_PI) * (111000);
                #else
                    f = 0;
                #endif
                
                pgH_px = model.gravity / dx_for_u * (model.h[p][i+1][j] - model.h[p][i-1][j]);

                pU2_px = 0.5 / dx_for_u * 
                        (model.gUpper[i+1][j][0] * pow(model..u[p][i+1][j], 2) - 
                         model.gUpper[i-1][j][0] * pow(model.u[p][i-1][j], 2));

                pUV_px = 1 / dx_for_u * 
                        (model.gUpper[i+1][j][1] * model..u[p][i+1][j] * model.v[p][i+1][j] - 
                         model.gUpper[i-1][j][1] * model.u[p][i-1][j] * model.v[p][i-1][j]);

                pV2_px = 0.5 / dx_for_u * 
                        (model.gUpper[i+1][j][3] * pow(model.v[p][i+1][j], 2) -
                         model.gUpper[i-1][j][3] * pow(model.v[p][i-1][j], 2));

                rotationU = (((model.v[p][i+1][j] - model.v[p][i-1][j]) / dx_for_u) - 
                             ((model.u[p][i][j+1] - model.u[p][i][j-1]) / dy_for_u) + model.sqrtG[i][j] * f) * 
                            (model.gUpper[i][j][2] * model.u[p][i][j] + model.gUpper[i][j][3] * model.v[p][i][j]);
            

                #ifdef Mountain
                    pgHs_px = model.gravity / dx_for_u * (model.hs[p][i+1][j] - model.csswm[p].hs[i-1][j]);
                    model.up[p][i][j] = model.um[p][i][j] + model.d2t * (-pgH_px - pgHs_px - pU2_px - pUV_px - pV2_px + rotationU);
                #else
                    model.up[p][i][j] = model.um[p][i][j] + model.d2t * (-pgH_px - pU2_px - pUV_px - pV2_px + rotationU);
                #endif

                #ifdef DIFFUSION
                    model.up[p][i][j] += model.d2t * KX * (model.um[p][i+1][j] - 2. * model.um[p][i][j] + model.um[p][i-1][j]) / pow(dx_for_u, 2) + 
                                               model.d2t * KY * (model.um[p][i][j+1] - 2. * model.um[p][i][j] + model.um[p][i][j-1]) / pow(dy_for_u, 2);
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
                dx_for_v = model.x[p][i+1][j] - model.x[p][i-1][j];
                dy_for_v = model.y[p][i][j+1] - model.y[p][i][j-1];

                #if defined(SteadyGeostrophy) || defined(Mountain)
                    f = 2 * OMEGA * (-cos(model.lon[p][i][j]) * cos(model.lat[p][i][j]) * sin(ALPHA0) + sin(model.lat[p][i][j]) * cos(ALPHA0));
                #elif defined(Barotropic) || defined(RossbyHaurwitz)
                    f = 2 * OMEGA * sin(model.lat[p][i][j]);
                #elif defined(EquatorialWave)
                    double f0 = 0;
                    double beta = 2.5 * 10E-5;
                    f = f0 + beta * abs(model.lat[p][i][j] * 180 / M_PI) * (111000);
                #else
                    f = 0;
                #endif

                pgH_py = model.gravity / dy_for_v * (model.h[p][i][j+1] - model.h[p][i][j-1]);

                pU2_py = 0.5 / dy_for_v * 
                        (model.gUpper[i][j+1][0] * pow(model.u[p][i][j+1], 2) - 
                         model.gUpper[i][j-1][0] * pow(model.u[p][i][j-1], 2));

                pUV_py = 1 / dy_for_v * 
                        (model.gUpper[i][j+1][1] * model.u[p][i][j+1] * model.v[p][i][j+1] - 
                         model.gUpper[i][j-1][1] * model.u[p][i][j-1] * model.v[p][i][j-1]);

                pV2_py = 0.5 / dy_for_v * 
                        (model.gUpper[i][j+1][3] * pow(model.v[p][i][j+1], 2) -
                         model.gUpper[i][j-1][3] * pow(model.v[p][i][j-1], 2));

                rotationV = (((model.v[p][i+1][j] - model.v[p][i-1][j]) / dx_for_v) - 
                             ((model.u[p][i][j+1] - model.u[p][i][j-1]) / dy_for_v) + model.sqrtG[i][j] * f) * 
                            (model.gUpper[i][j][0] * model.u[p][i][j] + model.gUpper[i][j][1] * model.v[p][i][j]);

                
                #ifdef Mountain
                    pgHs_py = model.gravity / dy_for_v * (model.csswm[p].hs[i][j+1] - model.csswm[p].hs[i][j-1]);;
                    model.vp[p][i][j] = model.vm[p][i][j] + model.d2t * (-pgH_py - pgHs_py - pU2_py - pUV_py - pV2_py - rotationV);
                #else
                    model.vp[p][i][j] = model.vm[p][i][j] + model.d2t * (-pgH_py - pU2_py - pUV_py - pV2_py - rotationV);
                #endif

                #ifdef DIFFUSION
                    model.vp[p][i][j] += model.d2t * KX * (model.vm[p][i+1][j] - 2. * model.vm[p][i][j] + model.vm[p][i-1][j]) / pow(dx_for_v, 2) + 
                                               model.d2t * KY * (model.vm[p][i][j+1] - 2. * model.vm[p][i][j] + model.vm[p][i][j-1]) / pow(dy_for_v, 2);
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
        for (int i = 2; i < model.nx-2; i++) {
            for (int j = 2; j < model.ny-2; j++) {
                dx_for_h = 0.5 * (model.x[p][i+1][j] - model.x[p][i-1][j]);
                dy_for_h = 0.5 * (model.y[p][i][j+1] - model.y[p][i][j-1]);

                psqrtGHU_px = (1. / (model.sqrtG[i][j] * 12. * dx_for_h)) * 
                              (-1.*(model.sqrtG[i+2][j] * model.h[p][i+2][j] * (model.gUpper[i+2][j][0] * model.u[p][i+2][j] + model.gUpper[i+2][j][1] * model.v[p][i+2][j]))
                               +8.*(model.sqrtG[i+1][j] * model.h[p][i+1][j] * (model.gUpper[i+1][j][0] * model.u[p][i+1][j] + model.gUpper[i+1][j][1] * model.v[p][i+1][j]))
                               -8.*(model.sqrtG[i-1][j] * model.h[p][i-1][j] * (model.gUpper[i-1][j][0] * model.u[p][i-1][j] + model.gUpper[i-1][j][1] * model.v[p][i-1][j]))
                               +1.*(model.sqrtG[i-2][j] * model.h[p][i-2][j] * (model.gUpper[i-2][j][0] * model.u[p][i-2][j] + model.gUpper[i-2][j][1] * model.v[p][i-2][j])));

                psqrtGHU_py = (1. / (model.sqrtG[i][j] * 12. * dy_for_h)) * 
                              (-1.*(model.sqrtG[i][j+2] * model.h[p][i][j+2] * (model.gUpper[i][j+2][2] * model.u[p][i][j+2] + model.gUpper[i][j+2][3] * model.v[p][i][j+2])) 
                               +8.*(model.sqrtG[i][j+1] * model.h[p][i][j+1] * (model.gUpper[i][j+1][2] * model.u[p][i][j+1] + model.gUpper[i][j+1][3] * model.v[p][i][j+1]))
                               -8.*(model.sqrtG[i][j-1] * model.h[p][i][j-1] * (model.gUpper[i][j-1][2] * model.u[p][i][j-1] + model.gUpper[i][j-1][3] * model.v[p][i][j-1]))
                               +1.*(model.sqrtG[i][j-2] * model.h[p][i][j-2] * (model.gUpper[i][j-2][2] * model.u[p][i][j-2] + model.gUpper[i][j-2][3] * model.v[p][i][j-2])));
            
                #if defined(AB2Time)
                    model.dh[p][i][j][(model.step+1)%2] = (-psqrtGHU_px - psqrtGHU_py);
                    #if defined(EquatorialWave)
                        if (model.status_add_forcing == true) model.dh[p][i][j][(model.step+1)%2] += model.h_forcing[p][i][j];
                    #endif
                    if (model.step == 0) model.dh[p][i][j][0] = model.dh[p][i][j][1];
                    model.hp[p][i][j] = model.h[p][i][j] + 1.5*model.dt*model.dh[p][i][j][(model.step+1)%2] - 0.5*model.dt*model.dh[p][i][j][model.step%2];
                #else
                    #if defined(EquatorialWave)
                        if (model.status_add_forcing == true) model.hp[p][i][j] = model.hm[p][i][j] + model.d2t * (-psqrtGHU_px - psqrtGHU_py + model.h_forcing[p][i][j]);
                        else model.hp[p][i][j] = model.hm[p][i][j] + model.d2t * (-psqrtGHU_px - psqrtGHU_py);
                    #else
                        model.hp[p][i][j] = model.hm[p][i][j] + model.d2t * (-psqrtGHU_px - psqrtGHU_py);
                    #endif
                #endif
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
        for (int i = 2; i < model.nx-2; i++) {
            for (int j = 2; j < model.ny-2; j++) {
                dx_for_u = 0.5 * (model.x[p][i+1][j] - model.x[p][i-1][j]);
                dy_for_u = 0.5 * (model.y[p][i][j+1] - model.y[p][i][j-1]);

                #if defined(SteadyGeostrophy) || defined(Mountain)
                    f = 2 * OMEGA * (-cos(model.lon[p][i][j]) * cos(model.lat[p][i][j]) * sin(ALPHA0) + sin(model.lat[p][i][j]) * cos(ALPHA0));
                #elif defined(Barotropic) || defined(RossbyHaurwitz)
                    f = 2 * OMEGA * sin(model.lat[p][i][j]);
                #elif defined(EquatorialWave) || defined(Uniform_f)
                    double f0 = 0;
                    double beta = 2.5 * 10E-11;
                    f = f0 + beta * model.lat[p][i][j] * 180. / M_PI * (111000.);
                #else
                    f = 0;
                #endif
                
                pgH_px = model.gravity / (12.*dx_for_u) * (-1*model.h[p][i+2][j] + 8*model.h[p][i+1][j] - 8*model.h[p][i-1][j] + 1*model.h[p][i-2][j]);

                pU2_px = 0.5 / (12.*dx_for_u) * 
                        (-1.*(model.gUpper[i+2][j][0] * pow(model.u[p][i+2][j], 2))
                         +8.*(model.gUpper[i+1][j][0] * pow(model.u[p][i+1][j], 2))
                         -8.*(model.gUpper[i-1][j][0] * pow(model.u[p][i-1][j], 2))
                         +1.*(model.gUpper[i-2][j][0] * pow(model.u[p][i-2][j], 2)));

                pUV_px = 1./ (12.*dx_for_u) * 
                        (-1.*(model.gUpper[i+2][j][1] * model.u[p][i+2][j] * model.v[p][i+2][j])
                         +8.*(model.gUpper[i+1][j][1] * model.u[p][i+1][j] * model.v[p][i+1][j])
                         -8.*(model.gUpper[i-1][j][1] * model.u[p][i-1][j] * model.v[p][i-1][j])
                         +1.*(model.gUpper[i-2][j][1] * model.u[p][i-2][j] * model.v[p][i-2][j]));

                pV2_px = 0.5 / (12.*dx_for_u) * 
                        (-1.*(model.gUpper[i+2][j][3] * pow(model.v[p][i+2][j], 2))
                         +8.*(model.gUpper[i+1][j][3] * pow(model.v[p][i+1][j], 2))
                         -8.*(model.gUpper[i-1][j][3] * pow(model.v[p][i-1][j], 2))
                         +1.*(model.gUpper[i-2][j][3] * pow(model.v[p][i-2][j], 2)));

                rotationU = (((-1.*model.v[p][i+2][j] + 8.*model.v[p][i+1][j] - 8.*model.v[p][i-1][j] + 1.*model.v[p][i-2][j]) / (12.*dx_for_u)) - 
                             ((-1.*model.u[p][i][j+2] + 8.*model.u[p][i][j+1] - 8.*model.u[p][i][j-1] + 1.*model.u[p][i][j-2]) / (12.*dy_for_u)) 
                             + model.sqrtG[i][j] * f) 
                             * (model.gUpper[i][j][2] * model.u[p][i][j] + model.gUpper[i][j][3] * model.v[p][i][j]);
            
                #if defined(AB2Time)
                    model.du[p][i][j][(model.step+1)%2] = (-pgH_px - pU2_px - pUV_px - pV2_px + rotationU);
                    #if defined(Mountain)
                        pgHs_px = model.gravity / (12.*dx_for_u) * (-1.*model.hs[p][i+2][j] + 8.*model.hs[p][i+1][j] - 8.*model.hs[p][i-1][j] + 1.*model.hs[p][i-2][j]);
                        model.du[p][i][j][(model.step+1)%2] += -pgHs_px;
                    #endif
                    if (model.step == 0) model.du[p][i][j][0] = model.du[p][i][j][1];
                    model.up[p][i][j] = model.u[p][i][j] + 1.5*model.dt*model.du[p][i][j][(model.step+1)%2] - 0.5*model.dt*model.du[p][i][j][model.step%2];
                #else
                    #if defined(Mountain)
                        pgHs_px = model.gravity / (12.*dx_for_u) * (-1.*model.hs[p][i+2][j] + 8.*model.hs[p][i+1][j] - 8.*model.hs[p][i-1][j] + 1.*model.hs[p][i-2][j]);
                        model.up[p][i][j] = model.um[p][i][j] + model.d2t * (-pgH_px - pgHs_px - pU2_px - pUV_px - pV2_px + rotationU);
                    #else
                        model.up[p][i][j] = model.um[p][i][j] + model.d2t * (-pgH_px - pU2_px - pUV_px - pV2_px + rotationU);
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
        for (int i = 2; i < model.nx-2; i++) {
            for (int j = 2; j < model.ny-2; j++) {
                dx_for_v = 0.5 * (model.x[p][i+1][j] - model.x[p][i-1][j]);
                dy_for_v = 0.5 * (model.y[p][i][j+1] - model.y[p][i][j-1]);

                #if defined(SteadyGeostrophy) || defined(Mountain)
                    f = 2 * OMEGA * (-cos(model.lon[p][i][j]) * cos(model.lat[p][i][j]) * sin(ALPHA0) + sin(model.lat[p][i][j]) * cos(ALPHA0));
                #elif defined(Barotropic) || defined(RossbyHaurwitz)
                    f = 2 * OMEGA * sin(model.lat[p][i][j]);
                #elif defined(EquatorialWave) || defined(Uniform_f)
                    double f0 = 0;
                    double beta = 2.5 * 10E-11;
                    f = f0 + beta * model.lat[p][i][j] * 180. / M_PI * (111000.);
                #else
                    f = 0;
                #endif

                pgH_py = model.gravity / (12.*dy_for_v) * (-1.*model.h[p][i][j+2] + 8.*model.h[p][i][j+1] - 8.*model.h[p][i][j-1] + 1.*model.h[p][i][j-2]);

                pU2_py = 0.5 /(12.*dy_for_v) * 
                        (-1.*(model.gUpper[i][j+2][0] * pow(model.u[p][i][j+2], 2))
                         +8.*(model.gUpper[i][j+1][0] * pow(model.u[p][i][j+1], 2))
                         -8.*(model.gUpper[i][j-1][0] * pow(model.u[p][i][j-1], 2))
                         +1.*(model.gUpper[i][j-2][0] * pow(model.u[p][i][j-2], 2)));

                pUV_py = 1. / (12.*dy_for_v) * 
                        (-1.*(model.gUpper[i][j+2][1] * model.u[p][i][j+2] * model.v[p][i][j+2])
                         +8.*(model.gUpper[i][j+1][1] * model.u[p][i][j+1] * model.v[p][i][j+1])
                         -8.*(model.gUpper[i][j-1][1] * model.u[p][i][j-1] * model.v[p][i][j-1])
                         +1.*(model.gUpper[i][j-2][1] * model.u[p][i][j-2] * model.v[p][i][j-2]));

                pV2_py = 0.5 / (12.*dy_for_v) * 
                        (-1.*(model.gUpper[i][j+2][3] * pow(model.v[p][i][j+2], 2))
                         +8.*(model.gUpper[i][j+1][3] * pow(model.v[p][i][j+1], 2))
                         -8.*(model.gUpper[i][j-1][3] * pow(model.v[p][i][j-1], 2))
                         +1.*(model.gUpper[i][j-2][3] * pow(model.v[p][i][j-2], 2)));

                rotationV = (((-1.*model.v[p][i+2][j] + 8.*model.v[p][i+1][j] - 8.*model.v[p][i-1][j] + 1.*model.v[p][i-2][j]) / (12.*dx_for_v)) - 
                             ((-1.*model.u[p][i][j+2] + 8.*model.u[p][i][j+1] - 8.*model.u[p][i][j-1] + 1.*model.u[p][i][j-2]) / (12.*dy_for_v)) 
                             + model.sqrtG[i][j] * f) * 
                            (model.gUpper[i][j][0] * model.u[p][i][j] + model.gUpper[i][j][1] * model.v[p][i][j]);

                #if defined(AB2Time)
                    model.dv[p][i][j][(model.step+1)%2] = (-pgH_py - pU2_py - pUV_py - pV2_py - rotationV);
                    #if defined(Mountain)
                        pgHs_py = model.gravity / (12.*dy_for_u) * (-1.*model.hs[p][i+2][j] + 8.*model.hs[p][i+1][j] - 8.*model.hs[p][i-1][j] + 1.*model.hs[p][i-2][j]);
                        model.dv[p][i][j][(model.step+1)%2] += -pgHs_py;
                    #endif
                    if (model.step == 0) model.dv[p][i][j][0] = model.dv[p][i][j][1];
                    model.vp[p][i][j] = model.v[p][i][j] + 1.5*model.dt*model.dv[p][i][j][(model.step+1)%2] - 0.5*model.dt*model.dv[p][i][j][model.step%2];
                #else
                    #if defined(Mountain)
                        pgHs_py = model.gravity / (12.*dy_for_v) * (-1.*model.hs[p][i][j+2] + 8.*model.hs[p][i][j+1] - 8.*model.hs[p][i][j-1] + 1.*model.hs[p][i][j-2]);
                        model.vp[p][i][j] = model.vm[p][i][j] + model.d2t * (-pgH_py - pgHs_py - pU2_py - pUV_py - pV2_py - rotationV);
                    #else
                        model.vp[p][i][j] = model.vm[p][i][j] + model.d2t * (-pgH_py - pU2_py - pUV_py - pV2_py - rotationV);
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
        for (int i = 0; i < model.nx; i++) {
            for (int j = 0; j < model.ny; j++) {
                model.hm[p][i][j] = model.h[p][i][j];
                model.h[p][i][j] = model.hp[p][i][j];

                model.um[p][i][j] = model.u[p][i][j];
                model.u[p][i][j] = model.up[p][i][j];

                model.vm[p][i][j] = model.v[p][i][j];
                model.v[p][i][j] = model.vp[p][i][j];
            }
        }
    }
    return;
}

void CSSWM::Iteration::TimeMarching(CSSWM &model) {
    CSSWM::Outputs::create_all_directory(model.outputpath);
    #ifdef NCOUTPUT
        CSSWM::Outputs::grid_nc(model);
    #endif
    #ifdef TXTOUTPUT
        CSSWM::Outputs::grid(model);
    #endif
    model.step = 0;
    // double timenow = 0.;
    double temp = model.timeend / model.dt;
    int nmax = (int) temp;

    while (model.step < nmax) {
        std::cout << model.step << std::endl;

        if (model.step % model.outputstep == 0) {
            #ifdef TXTOUTPUT
                CSSWM::Outputs::h(n, model);
                CSSWM::Outputs::u(n, model);
                CSSWM::Outputs::v(n, model);
            #endif

            #ifdef NCOUTPUT
                CSSWM::Outputs::huv_nc(model.step, model);
            #endif
        }

        #if defined(EquatorialWave)
            if (model.step * model.dt >= model.ADDFORCINGTIME) model.status_add_forcing = false;
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
        
        #if defined(TIMEFILTER) && !defined(AB2Time)
            CSSWM::NumericalProcess::timeFilterAll(model);
        #endif

        // next step
        CSSWM::Iteration::nextTimeStep(model);

        model.step++;
    }
    return;
}

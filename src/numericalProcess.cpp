#include "construction.hpp"

void CSSWM::NumericalProcess::DiffusionAll(CSSWM &model) {
    double dx, dy;
    for (int p = 0; p < 6; p++) {
        for (int i = 2; i < NX-2; i++) {
            for (int j = 2; j < NY-2; j++) {
                dx = 0.5 * (model.csswm[p].x[i+1][j] - model.csswm[p].x[i-1][j]);
                dy = 0.5 * (model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j-1]);

                model.csswm[p].hp[i][j] += model.d2t * model.diffusion_kx * (model.csswm[p].hm[i+1][j] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i-1][j]) / pow(dx, 2) + 
                                     model.d2t * model.diffusion_ky * (model.csswm[p].hm[i][j+1] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i][j-1]) / pow(dy, 2);
                model.csswm[p].up[i][j] += model.d2t * model.diffusion_kx * (model.csswm[p].um[i+1][j] - 2. * model.csswm[p].um[i][j] + model.csswm[p].um[i-1][j]) / pow(dx, 2) + 
                                     model.d2t * model.diffusion_ky * (model.csswm[p].um[i][j+1] - 2. * model.csswm[p].um[i][j] + model.csswm[p].um[i][j-1]) / pow(dy, 2);
                model.csswm[p].vp[i][j] += model.d2t * model.diffusion_kx * (model.csswm[p].vm[i+1][j] - 2. * model.csswm[p].vm[i][j] + model.csswm[p].vm[i-1][j]) / pow(dx, 2) + 
                                     model.d2t * model.diffusion_ky * (model.csswm[p].vm[i][j+1] - 2. * model.csswm[p].vm[i][j] + model.csswm[p].vm[i][j-1]) / pow(dy, 2);
            }
        }
    }
}

void CSSWM::NumericalProcess::timeFilterAll(CSSWM &model) {
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                model.csswm[p].h[i][j] += model.diffusion_ts * (model.csswm[p].hp[i][j] - 2 * model.csswm[p].h[i][j] + model.csswm[p].hm[i][j]);
                model.csswm[p].u[i][j] += model.diffusion_ts * (model.csswm[p].up[i][j] - 2 * model.csswm[p].u[i][j] + model.csswm[p].um[i][j]);
                model.csswm[p].v[i][j] += model.diffusion_ts * (model.csswm[p].vp[i][j] - 2 * model.csswm[p].v[i][j] + model.csswm[p].vm[i][j]);
            }
        }
    }
    return;
}

void CSSWM::NumericalProcess::NudgeH(CSSWM &model) {
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                model.csswm[p].hp[i][j] -= model.dt / (model.csswm_h_nudge_time) * (model.csswm[p].hp[i][j]-10454.608791605699);
            }
        }
    }
    return;
}
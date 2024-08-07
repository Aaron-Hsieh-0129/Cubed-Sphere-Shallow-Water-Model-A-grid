#include "construction.hpp"

void CSSWM::NumericalProcess::DiffusionAll(CSSWM &model) {
    double dx, dy;
    for (int p = 0; p < 6; p++) {
        for (int i = 2; i < model.nx-2; i++) {
            for (int j = 2; j < model.ny-2; j++) {
                dx = 0.5 * (model.x[p][i+1][j] - model.x[p][i-1][j]);
                dy = 0.5 * (model.y[p][i][j+1] - model.y[p][i][j-1]);

                model.hp[p][i][j] += model.d2t * model.Kx * (model.hm[p][i+1][j] - 2. * model.hm[p][i][j] + model.hm[p][i-1][j]) / pow(dx, 2) + 
                                     model.d2t * model.Ky * (model.hm[p][i][j+1] - 2. * model.hm[p][i][j] + model.hm[p][i][j-1]) / pow(dy, 2);
                model.up[p][i][j] += model.d2t * model.Kx * (model.um[p][i+1][j] - 2. * model.um[p][i][j] + model.um[p][i-1][j]) / pow(dx, 2) + 
                                     model.d2t * model.Ky * (model.um[p][i][j+1] - 2. * model.um[p][i][j] + model.um[p][i][j-1]) / pow(dy, 2);
                model.vp[p][i][j] += model.d2t * model.Kx * (model.vm[p][i+1][j] - 2. * model.vm[p][i][j] + model.vm[p][i-1][j]) / pow(dx, 2) + 
                                     model.d2t * model.Ky * (model.vm[p][i][j+1] - 2. * model.vm[p][i][j] + model.vm[p][i][j-1]) / pow(dy, 2);
            }
        }
    }
}

void CSSWM::NumericalProcess::timeFilterAll(CSSWM &model) {
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < model.nx; i++) {
            for (int j = 0; j < model.ny; j++) {
                model.h[p][i][j] += model.TIMETS * (model.hp[p][i][j] - 2 * model.h[p][i][j] + model.hm[p][i][j]);
                model.u[p][i][j] += model.TIMETS * (model.up[p][i][j] - 2 * model.u[p][i][j] + model.um[p][i][j]);
                model.v[p][i][j] += model.TIMETS * (model.vp[p][i][j] - 2 * model.v[p][i][j] + model.vm[p][i][j]);
            }
        }
    }
    return;
}

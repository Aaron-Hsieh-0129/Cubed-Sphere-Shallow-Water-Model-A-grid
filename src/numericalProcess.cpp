#include "construction.hpp"

void CSSWM::NumericalProcess::DiffusionAll(CSSWM &model) {
    double dx, dy;
    for (int p = 0; p < 6; p++) {
        for (int i = 2; i < NX-2; i++) {
            for (int j = 2; j < NY-2; j++) {
                dx = 0.5 * (model.x[p][i+1][j] - model.x[p][i-1][j]);
                dy = 0.5 * (model.y[p][i][j+1] - model.y[p][i][j-1]);

                model.hp[p][i][j] += D2T * KX * (model.hm[p][i+1][j] - 2. * model.hm[p][i][j] + model.hm[p][i-1][j]) / pow(dx, 2) + 
                                           D2T * KY * (model.hm[p][i][j+1] - 2. * model.hm[p][i][j] + model.hm[p][i][j-1]) / pow(dy, 2);
                model.up[p][i][j] += D2T * KX * (model.um[p][i+1][j] - 2. * model.um[p][i][j] + model.um[p][i-1][j]) / pow(dx, 2) + 
                                           D2T * KY * (model.um[p][i][j+1] - 2. * model.um[p][i][j] + model.um[p][i][j-1]) / pow(dy, 2);
                model.vp[p][i][j] += D2T * KX * (model.vm[p][i+1][j] - 2. * model.vm[p][i][j] + model.vm[p][i-1][j]) / pow(dx, 2) + 
                                           D2T * KY * (model.vm[p][i][j+1] - 2. * model.vm[p][i][j] + model.vm[p][i][j-1]) / pow(dy, 2);
            }
        }
    }
}

void CSSWM::NumericalProcess::timeFilterAll(CSSWM &model) {
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                model.h[p][i][j] += TIMETS * (model.hp[p][i][j] - 2 * model.h[p][i][j] + model.hm[p][i][j]);
                model.u[p][i][j] += TIMETS * (model.up[p][i][j] - 2 * model.u[p][i][j] + model.um[p][i][j]);
                model.v[p][i][j] += TIMETS * (model.vp[p][i][j] - 2 * model.v[p][i][j] + model.vm[p][i][j]);
            }
        }
    }
    return;
}

#include "construction.hpp"

void CSSWM::NumericalProcess::DiffusionAll(CSSWM &model) {
    double dx, dy;
    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < model.csswm[p].nx-1; i++) {
            for (int j = 1; j < model.csswm[p].ny-1; j++) {
                dx = 0.5 * (model.csswm[p].x[i+1][j] - model.csswm[p].x[i-1][j]);
                dy = 0.5 * (model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j-1]);

                model.csswm[p].hp[i][j] += model.d2t * model.KX * (model.csswm[p].hm[i+1][j] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i-1][j]) / pow(dx, 2) + 
                                           model.d2t * model.KY * (model.csswm[p].hm[i][j+1] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i][j-1]) / pow(dy, 2);
                model.csswm[p].up[i][j] += model.d2t * model.KX * (model.csswm[p].um[i+1][j] - 2. * model.csswm[p].um[i][j] + model.csswm[p].um[i-1][j]) / pow(dx, 2) + 
                                           model.d2t * model.KY * (model.csswm[p].um[i][j+1] - 2. * model.csswm[p].um[i][j] + model.csswm[p].um[i][j-1]) / pow(dy, 2);
                model.csswm[p].vp[i][j] += model.d2t * model.KX * (model.csswm[p].vm[i+1][j] - 2. * model.csswm[p].vm[i][j] + model.csswm[p].vm[i-1][j]) / pow(dx, 2) + 
                                           model.d2t * model.KY * (model.csswm[p].vm[i][j+1] - 2. * model.csswm[p].vm[i][j] + model.csswm[p].vm[i][j-1]) / pow(dy, 2);
            }
        }
    }
}

void CSSWM::NumericalProcess::timeFilterAll(CSSWM &model) {
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < model.nx; i++) {
            for (int j = 0; j < model.ny; j++) {
                model.csswm[p].h[i][j] += model.TIMETS * (model.csswm[p].hp[i][j] - 2 * model.csswm[p].h[i][j] + model.csswm[p].hm[i][j]);
                model.csswm[p].u[i][j] += model.TIMETS * (model.csswm[p].up[i][j] - 2 * model.csswm[p].u[i][j] + model.csswm[p].um[i][j]);
                model.csswm[p].v[i][j] += model.TIMETS * (model.csswm[p].vp[i][j] - 2 * model.csswm[p].v[i][j] + model.csswm[p].vm[i][j]);
            }
        }
    }
    return;
}

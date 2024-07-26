#include "construction.hpp"

int main(void) {
    clock_t start, stop;
    start = clock();

    Config_CSSWM config(30., 1., 1., 9.80665, 86400 * 3 * 24., "/data/Aaron/CSSWM_test/EquatorialWave/", 
                        1, 1E6, 1E6, 0.06, 1200. * 60.);
    CSSWM model(config);

    CSSWM::Init::Init2d(model);
    CSSWM::Iteration::TimeMarching(model);

    stop = clock();
    std::cout << double(stop - start) / CLOCKS_PER_SEC << " s" << std::endl;
    return 0;
}
#include "construction.hpp"

// Config_CSSWM config(180, 2, 2, 86400 * 3 * 24., 2E5, 2E5, 0.06, 9.80665, 10, 1200.*60.);
Config_CSSWM config(180, 2, 2, 86400 * 3 * 24., 2E5, 2E5, 0.06, 0.1, 10, 1200.*60., "/data/Aaron/TMIF/0527/csswm/");
CSSWM model(config);
int main(void) {
    clock_t start, stop;
    start = clock();

    CSSWM::Init::Init2d(model);
    CSSWM::Outputs::create_all_directory(model);
    CSSWM::Iteration::TimeMarching(model);

    stop = clock();
    std::cout << double(stop - start) / CLOCKS_PER_SEC << " s" << std::endl;
    return 0;
}
#include "iteration.hpp"

CSSWM model;
int main(void) {
    clock_t start, stop;
    start = clock();

    Init::Init2d(model);
    Iteration::leap_frog(model);

    stop = clock();
    std::cout << double(stop - start) / CLOCKS_PER_SEC << " s" << std::endl;
    return 0;
}
#include "construction.hpp"

CSSWM model;
int main(void) {
    clock_t start, stop;
    start = clock();

    CSSWM::Init::Init2d(model);
    CSSWM::Iteration::TimeMarching(model);

    stop = clock();
    std::cout << double(stop - start) / CLOCKS_PER_SEC << " s" << std::endl;
    return 0;
}
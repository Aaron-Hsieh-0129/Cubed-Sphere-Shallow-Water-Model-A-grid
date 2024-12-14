#include "construction.hpp"
#include "readConfig.hpp"

CSSWM model;
int main(void) {
    clock_t start, stop;
    start = clock();

    std::map<std::string, std::string> configs = csswm_read_config("../csswm_config.txt");
    std::string csswm_outputpath = configs["CSSWM_OUTPUTPATH"];
    double csswm_dt = std::stod(configs["CSSWM_DT"]);
    // double csswm_dx = std::stod(configs["CSSWM_DX"]);
    // double csswm_dy = std::stod(configs["CSSWM_DY"]);
    double csswm_timeend = std::stod(configs["CSSWM_TIMEEND"]);
    int csswm_outputstep = std::stoi(configs["CSSWM_OUTPUTSTEP"]);
    double csswm_gravity = std::stod(configs["CSSWM_GRAVITY"]);
    // int csswm_case = std::stoi(configs["CSSWM_CASE"]);
    double csswm_diffusion_kx = std::stod(configs["CSSWM_DIFFUSION_KX"]);
    double csswm_diffusion_ky = std::stod(configs["CSSWM_DIFFUSION_KY"]);
    double csswm_diffusion_ts = std::stod(configs["CSSWM_DIFFUSION_TS"]);

    // int csswm_nx = (int) (90 / csswm_dx + 2);
    // int csswm_ny = (int) (90 / csswm_dy + 2);
    double csswm_d2t = csswm_dt * csswm_dt;

    model.output_path = csswm_outputpath;
    model.dt = csswm_dt;
    model.d2t = csswm_d2t;
    // model.dx = csswm_dx;
    // model.dy = csswm_dy;
    // model.nx = csswm_nx;
    // model.ny = csswm_ny;
    model.timeend = csswm_timeend;
    model.outputstep = csswm_outputstep;
    model.gravity = csswm_gravity;
    // model.CASE = csswm_case;
    model.diffusion_kx = csswm_diffusion_kx;
    model.diffusion_ky = csswm_diffusion_ky;
    model.diffusion_ts = csswm_diffusion_ts;

    CSSWM::Init::Init2d(model);
    CSSWM::Iteration::TimeMarching(model);

    stop = clock();
    std::cout << double(stop - start) / CLOCKS_PER_SEC << " s" << std::endl;
    return 0;
}
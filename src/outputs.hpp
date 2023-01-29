#include <fstream>
#include <cstdlib>
#include <cstring>
#include <netcdf>
#include <vector>
#include "initialField.hpp"

class Outputs {
public:
    static void grid(CSSWM &);
    static void h(int, CSSWM &);
    static void u(int, CSSWM &);
    static void v(int, CSSWM &);

    static void grid_nc(CSSWM &);
    static void huv_nc(int, CSSWM &);

    static void create_all_directory();

private:
    static void create_directory(std::string);
};
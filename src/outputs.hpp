#include <fstream>
#include <cstdlib>
#include <cstring>
#include <netcdf>
#include <vector>
#include "initialField.hpp"

class Outputs {
public:
    static void output_parameter(CSSWM &);
    static void output_h(int, CSSWM &);
    static void output_u(int, CSSWM &);
    static void output_v(int, CSSWM &);

    static void output_parameter_nc(CSSWM &);

    static void create_all_directory();

private:
    static void create_directory(std::string);
};
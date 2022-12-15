#include <fstream>
#include <cstdlib>
#include <cstring>
#include "initialField.hpp"

class Outputs {
public:
    static void output_parameter(CSSWM &);
    static void output_h(int, CSSWM &);
    static void output_u(int, CSSWM &);
    static void output_v(int, CSSWM &);

private:
    static void create_directory(std::string);
};
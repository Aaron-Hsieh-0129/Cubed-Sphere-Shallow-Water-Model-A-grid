#include "constrcution.hpp"

class Init {
    public:
        Init();
        static void Init2d(CSSWM &);

    private:
        static double JungH(double, double);
        static double JungU(double, double);
        static double JungV(double);
        static double Gravity(double, double);
        static double SteadyGeostrophyH(double, double);
        static double ConvergenceRateH(double, double);
        static double BarotropicU(double);
        static double BarotropicH(double, double);
};
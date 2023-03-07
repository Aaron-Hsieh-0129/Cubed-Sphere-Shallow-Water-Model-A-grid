#include "construction.hpp"

class Init {
public:
    Init();
    static void Init2d(CSSWM &);
    static double BarotropicU(double);

private:
    static double JungH(double, double);
    static double JungU(double, double);
    static double JungV(double);
    
    static double AdvectionH(double, double);
    static double AdvectionU(double, double);
    static double AdvectionV(double);

    static double Gravity(double, double);
    static double SteadyGeostrophyH(double, double);
    static double SteadyGeostrophyU(double, double);
    static double SteadyGeostrophyV(double);
    static double ConvergenceRateH(double, double);
    static double DeformationalFlowH(double, double);
    static double BarotropicH(double);
    static double BarotropicHPrime(double, double);

    static double MountainH(double, double);
    static double MountainU(double, double);
    static double MountainV(double);

    static double RossbyHaurwitzH(double, double);
    static double RossbyHaurwitzU(double, double);
    static double RossbyHaurwitzV(double, double);

    static double simpson(double, double);
    static double func(double);
};
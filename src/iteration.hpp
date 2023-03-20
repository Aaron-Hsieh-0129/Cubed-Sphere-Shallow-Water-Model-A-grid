#include "outputs.hpp"

class Iteration {
public:
    
    static void leap_frog(CSSWM &);
    #if defined(SecondOrderSpace)
        static void ph_pt_2(CSSWM &);
        static void pu_pt_2(CSSWM &);
        static void pv_pt_2(CSSWM &);
    #elif defined(FourthOrderSpace) 
        static void ph_pt_4(CSSWM &);
        static void pu_pt_4(CSSWM &);
        static void pv_pt_4(CSSWM &);
    #endif
};
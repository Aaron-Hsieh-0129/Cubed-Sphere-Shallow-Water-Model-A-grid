#ifndef __DEFINE_H__
#define __DEFINE_H__

#define FILLVALUE (-99999999999999);

#define RADIUS (6371220.)
#define GRAVITY (9.80665)
#define OMEGA (7.292E-5)
#define BETA (2.5E-11)

#define TIMEEND (86400 * 3 * 24)
#define OUTPUTPATH "/data/Aaron/CSSWM_test/EquatorialWave/"
#define OUTPUTINTERVAL (1)
// #define SecondOrderSpace
#define FourthOrderSpace
#define NCOUTPUT
// #define TXTOUTPUT


#define D2T (2. * DT)

// Jung
#define ALPHA0 (0)
// #define Advection
// #define GravityWave
// #define SteadyGeostrophy
// #define Barotropic
// #define Mountain
// #define RossbyHaurwitz
#define EquatorialWave
// #define Uniform
// #define Uniform_f

// #define TrueSol
#define DIFFUSION
#define TIMEFILTER
#define KX (1E6)
#define KY (1E6)
#define TIMETS (0.06)

#define ADDFORCINGTIME (1200. * 60.)

#endif

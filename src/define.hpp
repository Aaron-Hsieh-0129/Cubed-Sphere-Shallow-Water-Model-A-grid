#ifndef __DEFINE_H__
#define __DEFINE_H__

#define FILLVALUE (-99999999999999);

#define RADIUS (6371220.)
#define GRAVITY (9.80665)
#define OMEGA (7.292E-5)
#define BETA (2.5E-11)

#define DX (1)
#define DY (1)
#define DT (200.)
#define TIMEEND (60000.)
#define OUTPUTPATH "/data/mog/TMIF_CSSWM/mog_test/csswm/"
#define OUTPUTINTERVAL (1)
// #define SecondOrderSpace
#define FourthOrderSpace
#define AB2Time
#define NCOUTPUT
// #define TXTOUTPUT

#if defined(SecondOrderSpace)
    #define NX ((int) (90/DX + 2))
    #define NY ((int) (90/DY + 2))
#elif defined(FourthOrderSpace) 
    #define NX ((int) (90/DX + 4))
    #define NY ((int) (90/DY + 4))
#endif

#define D2T (2. * DT)

#define ALPHA0 (0) // This is the angle of the wind direction
// Jung
// #define Advection
// #define GravityWave
// #define SteadyGeostrophy
// #define Barotropic
// #define Mountain
// #define RossbyHaurwitz // In this case, the gravity wave speed should be 9.80665
// #define EquatorialWave
#define Uniform

// #define TrueSol
#define DIFFUSION
#define TIMEFILTER
#define KX (2E5)
#define KY (2E5)
#define TIMETS (0.06)

#define ADDFORCINGTIME (1200. * 60.)

#endif

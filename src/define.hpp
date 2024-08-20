#ifndef __DEFINE_H__
#define __DEFINE_H__

#define FILLVALUE (-99999999999999);

#define RADIUS (6371220.)
#define GRAVITY (0.1)
#define OMEGA (7.292E-5)
#define BETA (2.5E-11)

#define DX (1)
#define DY (1)
#define DT (900.)
#define TIMEEND (86400 * 3 * 24)
#define OUTPUTPATH "/data/Aaron/TMIF/SP_old_CSSWM/AB2_prof_RKM_dt900_3_csswm_3_vvm_2E5diff_lineless_3B/csswm/"
// #define OUTPUTPATH "/data/Aaron/TMIF_CSSWM/sine2/csswm/"
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

// Jung
#define ALPHA0 (0)
// #define Advection
// #define GravityWave
// #define SteadyGeostrophy
// #define Barotropic
// #define Mountain
// #define RossbyHaurwitz
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

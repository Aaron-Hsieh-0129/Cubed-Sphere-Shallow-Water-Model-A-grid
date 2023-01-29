#ifndef __DEFINE_H__
#define __DEFINE_H__

#define FILLVALUE (-99999999999999);

#define RADIUS (6371220.)
#define GRAVITY (9.80616)
#define OMEGA (7.292E-5)

#define DX (2)
#define DY (2)
#define NX ((int) (90/DX + 2))
#define NY ((int) (90/DY + 2))
#define DT (180.)
#define D2T (2. * DT)
#define TIMEEND (86400 * 12 * 2)
#define OUTPUTINTERVAL (50)
#define OUTPUTPATH "../outputs/"
#define TXTOUTPUT
#define NCOUTPUT

#define ALPHA0 (0)
// #define ADVECTION
// #define Jung
// #define GravityWave
// #define SteadyGeostrophy
// #define Barotropic
// #define Mountain
#define RossbyHaurwitz

#define DIFFUSION
#define TIMEFILTER
#define KX (200000.)
#define KY (200000.)
#define TIMETS (0.06)

#endif
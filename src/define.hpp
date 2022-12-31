#ifndef __DEFINE_H__
#define __DEFINE_H__

#define FILLVALUE (-99999999999999);

#define RADIUS (6371220.)
#define GRAVITY (9.80616)
#define OMEGA (7.292E-5)

#define DX (0.5)
#define DY (0.5)
#define NX ((int) (90/DX + 2))
#define NY ((int) (90/DY + 2))
#define DT (45.)
#define D2T (2. * DT)
#define TIMEEND (86400 * 12 * 2)
#define OUTPUTINTERVAL (50)
#define OUTPUTPATH "../outputs/"

#define ALPHA0 (M_PI/4.)
// #define ADVECTION
// #define Jung
// #define GravityWave
// #define SteadyGeostrophy
#define Barotropic

#define DIFFUSION
#define TIMEFILTER
#define KX (150000.)
#define KY (150000.)
#define TIMETS (0.06)

#endif
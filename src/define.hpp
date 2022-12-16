#ifndef __DEFINE_H__
#define __DEFINE_H__

#define FILLVALUE (-99999999999999);

#define RADIUS (6371220.)
#define GRAVITY (9.80616)
#define OMEGA (7.292E-5)
#define CF (0.)

#define DX (2)
#define DY (2)
#define NX ((int) 90/DX + 2)
#define NY ((int) 90/DY + 2)
#define DT (360.)
#define D2T (2. * DT)
#define TIMEEND (24. * 86400 * 6)
#define OUTPUTINTERVAL (50)

#define ALPHA0 (M_PI/4.)
// #define Jung
#define GravityWave
// #define SteadyGeostrophy

#define DIFFUSION
#define TIMEFILTER
#define KX (1000000.)
#define KY (1000000.)
#define TIMETS (0.06)

#endif
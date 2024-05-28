#include <cmath>
#include <iostream>
#include <iomanip>
#include <memory>
#include "define.hpp"

class Config_CSSWM {
public:
    Config_CSSWM(double dt, double dx, double dy, double TIMEEND, double KX, double KY, double TIMETS, double GRAVITY, int OUTPUTSTEP, double ADDFORCINGTIME, std::string OUTPUTPATH) 
      : dt(dt), dx(dx), dy(dy), TIMEEND(TIMEEND), KX(KX), KY(KY), TIMETS(TIMETS), GRAVITY(GRAVITY), OUTPUTSTEP(OUTPUTSTEP), ADDFORCINGTIME(ADDFORCINGTIME), OUTPUTPATH(OUTPUTPATH) {}
    ~Config_CSSWM() = default;

    double dt;
    double dx;
    double dy;
    double TIMEEND;
    double KX;
    double KY;
    double TIMETS;
    double GRAVITY;
    int OUTPUTSTEP;
    double ADDFORCINGTIME;
    std::string OUTPUTPATH;
};

class CSSWM {
public:
    class patch {
    public:
        patch(const Config_CSSWM &config) 
          :  
            #if defined(SecondOrderSpace)
                nx(90/config.dx+2), ny(90/config.dx+2) 
            #elif defined(FourthOrderSpace)
                nx(90/config.dx+4), ny(90/config.dx+4) 
            #endif
        {
            allocateMemory();
            initializeArrays();
        }
        ~patch() {
            printf("Free patch\n");
            deallocateMemory();
        }

        void allocateMemory() {
            hp = allocate2DArray(nx, ny);
            h = allocate2DArray(nx, ny);
            hm = allocate2DArray(nx, ny);
            up = allocate2DArray(nx, ny);
            u = allocate2DArray(nx, ny);
            um = allocate2DArray(nx, ny);
            vp = allocate2DArray(nx, ny);
            v = allocate2DArray(nx, ny);
            vm = allocate2DArray(nx, ny);
            lon = allocate2DArray(nx, ny);
            lat = allocate2DArray(nx, ny);
            lon_original = allocate2DArray(nx, ny);
            x = allocate2DArray(nx, ny);
            y = allocate2DArray(nx, ny);
            A = allocate3DArray(nx, ny, 4);
            IA = allocate3DArray(nx, ny, 4);

            hpcont = new double[nx * ny];
            hcont = new double[nx * ny];
            hmcont = new double[nx * ny];
            upcont = new double[nx * ny];
            ucont = new double[nx * ny];
            umcont = new double[nx * ny];
            vpcont = new double[nx * ny];
            vcont = new double[nx * ny];
            vmcont = new double[nx * ny];

            loncont = new double[nx * ny];
            latcont = new double[nx * ny];
            lon_originalcont = new double[nx * ny];
            xcont = new double[nx * ny];
            ycont = new double[nx * ny];
            Acont = new double[nx * ny * 4];
            IAcont = new double[nx * ny * 4];

            #if defined(Mountain)
                hs = allocate2DArray(nx, ny);
                hscont = new double[nx * ny];
            #endif

            #if defined(EquatorialWave)
                h_forcing = allocate2DArray(nx, ny);
                h_forcingcont = new double[nx * ny];
            #endif

            #if defined(SecondOrderSpace)
                IP1_L = allocate2DArray(nx, 4);
                IP1_R = allocate2DArray(nx, 4);
                IP1_U = allocate2DArray(nx, 4);
                IP1_D = allocate2DArray(nx, 4);

                IP1_Lcont = new double[nx * 4];
                IP1_Rcont = new double[nx * 4];
                IP1_Ucont = new double[nx * 4];
                IP1_Dcont = new double[nx * 4];
            #elif defined(FourthOrderSpace)
                IP_ouTTer_L = allocate2DArray(nx, 4);
                IP_ouTTer_R = allocate2DArray(nx, 4);
                IP_ouTTer_U = allocate2DArray(nx, 4);
                IP_ouTTer_D = allocate2DArray(nx, 4);
                IP_ouTer_L = allocate2DArray(nx, 4);
                IP_ouTer_R = allocate2DArray(nx, 4);
                IP_ouTer_U = allocate2DArray(nx, 4);
                IP_ouTer_D = allocate2DArray(nx, 4);

                IP_ouTTer_Lcont = new double[nx * 4];
                IP_ouTTer_Rcont = new double[nx * 4];
                IP_ouTTer_Ucont = new double[nx * 4];
                IP_ouTTer_Dcont = new double[nx * 4];
                IP_ouTer_Lcont = new double[nx * 4];
                IP_ouTer_Rcont = new double[nx * 4];
                IP_ouTer_Ucont = new double[nx * 4];
                IP_ouTer_Dcont = new double[nx * 4];
            #endif
        }

        void initializeArrays() {
            for (int i = 0; i < nx; ++i) {
                hp[i] = &hpcont[i * ny];
                h[i] = &hcont[i * ny];
                hm[i] = &hmcont[i * ny];
                up[i] = &upcont[i * ny];
                u[i] = &ucont[i * ny];
                um[i] = &umcont[i * ny];
                vp[i] = &vpcont[i * ny];
                v[i] = &vcont[i * ny];
                vm[i] = &vmcont[i * ny];
                #if defined(Mountain)
                    hs[i] = &hscont[i * ny];
                #endif

                lon[i] = &loncont[i * ny];
                lat[i] = &latcont[i * ny];
                lon_original[i] = &lon_originalcont[i * ny];
                x[i] = &xcont[i * ny];
                y[i] = &ycont[i * ny];

                A[i] = allocate2DArray(ny, 4);
                IA[i] = allocate2DArray(ny, 4);
                for (int j = 0; j < ny; ++j) {
                    A[i][j] = &Acont[(i * ny + j) * 4];
                    IA[i][j] = &IAcont[(i * ny + j) * 4];
                }

                #if defined(SecondOrderSpace)
                    IP1_L[i] = &IP1_Lcont[i * 4];
                    IP1_R[i] = &IP1_Rcont[i * 4];
                    IP1_U[i] = &IP1_Ucont[i * 4];
                    IP1_D[i] = &IP1_Dcont[i * 4];
                #elif defined(FourthOrderSpace)
                    IP_ouTTer_L[i] = &IP_ouTTer_Lcont[i * 4];
                    IP_ouTTer_R[i] = &IP_ouTTer_Rcont[i * 4];
                    IP_ouTTer_U[i] = &IP_ouTTer_Ucont[i * 4];
                    IP_ouTTer_D[i] = &IP_ouTTer_Dcont[i * 4];
                    IP_ouTer_L[i] = &IP_ouTer_Lcont[i * 4];
                    IP_ouTer_R[i] = &IP_ouTer_Rcont[i * 4];
                    IP_ouTer_U[i] = &IP_ouTer_Ucont[i * 4];
                    IP_ouTer_D[i] = &IP_ouTer_Dcont[i * 4];
                #endif
            }

            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    hp[i][j] = h[i][j] = hm[i][j] = FILLVALUE;
                    up[i][j] = u[i][j] = um[i][j] = FILLVALUE;
                    vp[i][j] = v[i][j] = vm[i][j] = FILLVALUE;
                    lon[i][j] = lat[i][j] = lon_original[i][j] = FILLVALUE;
                    x[i][j] = y[i][j] = FILLVALUE;
                    for (int k = 0; k < 4; ++k) {
                        A[i][j][k] = IA[i][j][k] = FILLVALUE;
                    }
                    #if defined(Mountain)
                        hs[i][j] = FILLVALUE;
                    #endif
                }
            }
        }

        void deallocateMemory() {
            // for (int i = 0; i < nx; ++i) {
            //     delete[] hp[i];
            //     delete[] h[i];
            //     delete[] hm[i];
            //     delete[] up[i];
            //     delete[] u[i];
            //     delete[] um[i];
            //     delete[] vp[i];
            //     delete[] v[i];
            //     delete[] vm[i];
            //     delete[] lon[i];
            //     delete[] lat[i];
            //     delete[] lon_original[i];
            //     delete[] x[i];
            //     delete[] y[i];
            //     for (int j = 0; j < ny; ++j) {
            //         delete[] A[i][j];
            //         delete[] IA[i][j];
            //     }
            //     delete[] A[i];
            //     delete[] IA[i];

            delete[] hp;
            delete[] h;
            delete[] hm;
            delete[] up;
            delete[] u;
            delete[] um;
            delete[] vp;
            delete[] v;
            delete[] vm;
            delete[] lon;
            delete[] lat;
            delete[] lon_original;
            delete[] x;
            delete[] y;
            delete[] A;
            delete[] IA;

            delete[] hpcont;
            delete[] hcont;
            delete[] hmcont;
            delete[] upcont;
            delete[] ucont;
            delete[] umcont;
            delete[] vpcont;
            delete[] vcont;
            delete[] vmcont;
            delete[] loncont;
            delete[] latcont;
            delete[] lon_originalcont;
            delete[] xcont;
            delete[] ycont;
            delete[] Acont;
            delete[] IAcont;

            #if defined(Mountain)
                delete[] hscont;
            #endif

            #if defined(EquatorialWave)
                delete[] h_forcingcont;
            #endif

            #if defined(SecondOrderSpace)
                delete[] IP1_L;
                delete[] IP1_R;
                delete[] IP1_U;
                delete[] IP1_D;
                delete[] IP1_Lcont;
                delete[] IP1_Rcont;
                delete[] IP1_Ucont;
                delete[] IP1_Dcont;
            #elif defined(FourthOrderSpace)
                delete[] IP_ouTTer_L;
                delete[] IP_ouTTer_R;
                delete[] IP_ouTTer_U;
                delete[] IP_ouTTer_D;
                delete[] IP_ouTer_L;
                delete[] IP_ouTer_R;
                delete[] IP_ouTer_U;
                delete[] IP_ouTer_D;
                delete[] IP_ouTTer_Lcont;
                delete[] IP_ouTTer_Rcont;
                delete[] IP_ouTTer_Ucont;
                delete[] IP_ouTTer_Dcont;
                delete[] IP_ouTer_Lcont;
                delete[] IP_ouTer_Rcont;
                delete[] IP_ouTer_Ucont;
                delete[] IP_ouTer_Dcont;
            #endif
        }

        double** allocate2DArray(int dim1, int dim2) {
            double** array = new double*[dim1];
            for (int i = 0; i < dim1; ++i) {
                array[i] = new double[dim2];
            }
            return array;
        }

        double*** allocate3DArray(int dim1, int dim2, int dim3) {
            double*** array = new double**[dim1];
            for (int i = 0; i < dim1; ++i) {
                array[i] = new double*[dim2];
                for (int j = 0; j < dim2; ++j) {
                    array[i][j] = new double[dim3];
                }
            }
            return array;
        }

        int nx;
        int ny;

        double **hp, **h, **hm, **up, **u, **um, **vp, **v, **vm;
        double **lon, **lat, **lon_original, **x, **y;
        double ***A, ***IA;

        double *hpcont, *hcont, *hmcont, *upcont, *ucont, *umcont, *vpcont, *vcont, *vmcont;
        double *loncont, *latcont, *lon_originalcont, *xcont, *ycont, *Acont, *IAcont;

        #if defined(Mountain)
            double **hs;
            double *hscont;
        #endif

        #if defined(EquatorialWave)
            double **h_forcing;
            double *h_forcingcont;
        #endif

        #if defined(SecondOrderSpace)
            double **IP1_L, **IP1_R, **IP1_U, **IP1_D;
            double *IP1_Lcont, *IP1_Rcont, *IP1_Ucont, *IP1_Dcont;
        #elif defined(FourthOrderSpace)
            double **IP_ouTTer_L, **IP_ouTTer_R, **IP_ouTTer_U, **IP_ouTTer_D;
            double **IP_ouTer_L, **IP_ouTer_R, **IP_ouTer_U, **IP_ouTer_D;
            double *IP_ouTTer_Lcont, *IP_ouTTer_Rcont, *IP_ouTTer_Ucont, *IP_ouTTer_Dcont;
            double *IP_ouTer_Lcont, *IP_ouTer_Rcont, *IP_ouTer_Ucont, *IP_ouTer_Dcont;
        #endif
    };
    CSSWM(const Config_CSSWM &config) 
    : dt(config.dt), dx(config.dx), dy(config.dy), 
        #if defined(SecondOrderSpace)
            nx(90/config.dx+2), ny(90/config.dx+2),
        #elif defined(FourthOrderSpace)
            nx(90/config.dx+4), ny(90/config.dx+4),
        #endif
      TIMEEND(config.TIMEEND),
      KX(config.KX), KY(config.KY), TIMETS(config.TIMETS),
      GRAVITY(config.GRAVITY),
      OUTPUTSTEP(config.OUTPUTSTEP), ADDFORCINGTIME(config.ADDFORCINGTIME), OUTPUTPATH(config.OUTPUTPATH),
      d2t(2. * config.dt), 
      csswm{patch(config), patch(config), patch(config), patch(config), patch(config), patch(config)} {
        allocateMemory();
        initializeArrays();
        assignMatrix();
    }
    ~CSSWM() {
        printf("Free CSSWM\n");
        deallocateMemory();
    }

    void allocateMemory() {
        sqrtG = allocate2DArray(nx, ny);
        gamma = allocate2DArray(nx, ny);
        gLower = allocate3DArray(nx, ny, 4);
        gUpper = allocate3DArray(nx, ny, 4);
        alpha2D = allocate2DArray(nx, ny);
        beta2D = allocate2DArray(nx, ny);

        sqrtGcont = new double[nx * ny];
        gammacont = new double[nx * ny];
        gLowercont = new double[nx * ny * 4];
        gUppercont = new double[nx * ny * 4];
        alpha2Dcont = new double[nx * ny];
        beta2Dcont = new double[nx * ny];

        #if defined(SecondOrderSpace)
            match = allocate2DArrayInt(24, 8);
            checkIP = allocate2DArrayInt(nx, 2);

            matchcont = new int[24 * 8];
            checkIPcont = new int[nx * 2];
        #elif defined(FourthOrderSpace)
            match_ouTTer = allocate2DArrayInt(24, 8);
            match_ouTer = allocate2DArrayInt(24, 8);
            checkIP_ouTTer = allocate2DArrayInt(nx, 2);
            checkIP_ouTer = allocate2DArrayInt(nx, 2);

            match_ouTTercont = new int[24 * 8];
            match_ouTercont = new int[24 * 8];
            checkIP_ouTTercont = new int[nx * 2];
            checkIP_ouTercont = new int[nx * 2];
        #endif
    }

    void initializeArrays() {
        for (int i = 0; i < nx; i++) {
            sqrtG[i] = &sqrtGcont[i * ny];
            gamma[i] = &gammacont[i * ny];
            alpha2D[i] = &alpha2Dcont[i * ny];
            beta2D[i] = &beta2Dcont[i * ny];

            gLower[i] = allocate2DArray(ny, 4);
            gUpper[i] = allocate2DArray(ny, 4);
            for (int j = 0; j < ny; ++j) {
                gLower[i][j] = &gLowercont[(i * ny + j) * 4];
                gUpper[i][j] = &gUppercont[(i * ny + j) * 4];
            }
        }

        for (int i = 0; i < 24; i++) {
            #if defined(SecondOrderSpace)
                match[i] = &matchcont[i * 8];
                checkIP[i] = &checkIPcont[i * 2];
            #elif defined(FourthOrderSpace)
                match_ouTTer[i] = &match_ouTTercont[i * 8];
                match_ouTer[i] = &match_ouTercont[i * 8];
            #endif
        }
    }

    void assignMatrix() {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                #if defined(SecondOrderSpace)
                    alpha2D[i][j] = -M_PI/4. + (M_PI/2.) / (nx-2.) * (i-0.5);
                #elif defined(FourthOrderSpace)
                    alpha2D[i][j] = -M_PI/4. + (M_PI/2.) / (nx-4.) * (i-1.5);
                #endif
            }
        }
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                #if defined(SecondOrderSpace)
                    beta2D[i][j] = -M_PI/4. + (M_PI/2.) / (nx-2.) * (j-0.5);
                #elif defined(FourthOrderSpace)
                    beta2D[i][j] = -M_PI/4. + (M_PI/2.) / (nx-4.) * (j-1.5);
                #endif
            }
        }

        Construct_gamma_sqrtG_GUpper(alpha2D, beta2D, gamma, sqrtG, gUpper, gLower);

        for (int p = 0; p < 6; p++) {
            if (p == 4) {
                Construct_p4_lonlat_xy_AIA(alpha2D, beta2D, gamma, csswm[p].lon, csswm[p].lat, csswm[p].lon_original, csswm[p].x, csswm[p].y, csswm[p].A, csswm[p].IA);
                continue;
            }
            if (p == 5) {
                Construct_p5_lonlat_xy_AIA(alpha2D, beta2D, gamma, csswm[p].lon, csswm[p].lat, csswm[p].lon_original, csswm[p].x, csswm[p].y, csswm[p].A, csswm[p].IA);
                continue;
            }
            Construct_p0123_lonlat_xy_AIA(p, alpha2D, beta2D, gamma, csswm[p].lon, csswm[p].lat, csswm[p].lon_original, csswm[p].x, csswm[p].y, csswm[p].A, csswm[p].IA);
        }

        #if defined(SecondOrderSpace)
            // Interpolation matrix construction
            int idx = 0, ipIdx = 0;
            double A1, A2, B;

            // Construct Interpolation 2D array filled with index (Note that all interpolation relation between every boundary is same)
            while (idx < nx && ipIdx < nx - 1) {
                B = csswm[0].lat[nx-1][idx];
                A1 = csswm[1].lat[1][ipIdx], A2 = csswm[1].lat[1][ipIdx+1];

                if (A1 < B && B < A2) {
                    checkIP[idx][0] = ipIdx;
                    checkIP[idx][1] = ipIdx + 1;
                    idx++;
                }
                else if (A1 == B) {
                    checkIP[idx][0] = ipIdx;
                    checkIP[idx][1] = ipIdx;
                    idx++;
                }
                else if (A2 == B) {
                    checkIP[idx][0] = ipIdx + 1;
                    checkIP[idx][1] = ipIdx + 1;
                    idx++;
                }
                else {
                    ipIdx++;
                }
            } 
        #elif defined(FourthOrderSpace)
            // For ouTTer
            // Interpolation matrix construction
            int idx = 0, ipIdx = 0;
            double A1, A2, B;
            // Construct Interpolation 2D array filled with index (Note that all interpolation relation between every boundary is same)
            while (idx < nx && ipIdx < nx - 1) {
                B = csswm[0].lat[nx-1][idx];
                A1 = csswm[1].lat[3][ipIdx], A2 = csswm[1].lat[3][ipIdx+1];

                if (A1 < B && B < A2) {
                    checkIP_ouTTer[idx][0] = ipIdx;
                    checkIP_ouTTer[idx][1] = ipIdx + 1;
                    idx++;
                }
                else if (A1 == B) {
                    checkIP_ouTTer[idx][0] = ipIdx;
                    checkIP_ouTTer[idx][1] = ipIdx;
                    idx++;
                }
                else if (A2 == B) {
                    checkIP_ouTTer[idx][0] = ipIdx + 1;
                    checkIP_ouTTer[idx][1] = ipIdx + 1;
                    idx++;
                }
                else {
                    ipIdx++;
                }
            } 

            // For ouTer
            // Interpolation matrix construction
            idx = 0, ipIdx = 0;
            A1 = 0, A2 = 0, B = 0;
            // Construct Interpolation 2D array filled with index (Note that all interpolation relation between every boundary is same)
            while (idx < nx && ipIdx < nx - 1) {
                B = csswm[0].lat[nx-2][idx];
                A1 = csswm[1].lat[2][ipIdx], A2 = csswm[1].lat[2][ipIdx+1];

                if (A1 < B && B < A2) {
                    checkIP_ouTer[idx][0] = ipIdx;
                    checkIP_ouTer[idx][1] = ipIdx + 1;
                    idx++;
                }
                else if (A1 == B) {
                    checkIP_ouTer[idx][0] = ipIdx;
                    checkIP_ouTer[idx][1] = ipIdx;
                    idx++;
                }
                else if (A2 == B) {
                    checkIP_ouTer[idx][0] = ipIdx + 1;
                    checkIP_ouTer[idx][1] = ipIdx + 1;
                    idx++;
                }
                else {
                    ipIdx++;
                }
            }
        #endif

        // Construct a array for dealing with interpolation between all patch
        // (p1, p2, i1, j1, i2, j2, reversed, LonLat: 1/0) 
        // Left, Right, Up, Down
        #if defined(SecondOrderSpace)
            initMatch_1point(match);
        #elif defined(FourthOrderSpace)
            initMatch_2point(match_ouTer, match_ouTTer);
        #endif

        // Construct patch to patch transformation matrix (Declare at transform.cpp)
        Cube2Cube_matrix();

        // define mountain
        #ifdef Mountain
            for (int p = 0; p < 6; p++) {
                for (int i = 0; i < nx; i++) {
                    for (int j = 0; j < ny; j++) {
                        double r0 = M_PI / 9.;
                        double lonC = 3. * M_PI / 2., latC = M_PI / 6.;
                        double where = sqrt(pow(csswm[p].lon[i][j] - lonC, 2) + pow(csswm[p].lat[i][j] - latC, 2));
                        double r = r0 >  where ? where : r0;
                        double hs0 = 2000.;
                        csswm[p].hs[i][j] = hs0 * (1 - r / r0);
                    }
                }
            }
        #endif
    }

    void deallocateMemory() {
        // for (int i = 0; i < nx; ++i) {
        //     delete[] sqrtG[i];
        //     delete[] gamma[i];
        //     delete[] alpha2D[i];
        //     delete[] beta2D[i];
        //     for (int j = 0; j < ny; ++j) {
        //         delete[] gLower[i][j];
        //         delete[] gUpper[i][j];
        //     }
        //     delete[] gLower[i];
        //     delete[] gUpper[i];

        //     #if defined(SecondOrderSpace)
        //         delete[] checkIP[i];
        //     #elif defined(FourthOrderSpace)
        //         delete[] checkIP_ouTTer[i];
        //         delete[] checkIP_ouTer[i];
        //     #endif
        // }
        delete[] sqrtG;
        delete[] gamma;
        delete[] alpha2D;
        delete[] beta2D;
        delete[] gLower;
        delete[] gUpper;

        // for (int i = 0; i < 24; ++i) {
        //     #if defined(SecondOrderSpace)
        //         delete[] match[i];
        //     #elif defined(FourthOrderSpace)
        //         delete[] match_ouTTer[i];
        //         delete[] match_ouTer[i];
        //     #endif
        // }

        delete[] sqrtGcont;
        delete[] gammacont;
        delete[] gLowercont;
        delete[] gUppercont;
        delete[] alpha2Dcont;
        delete[] beta2Dcont;

        #if defined(SecondOrderSpace)
            delete[] match;
            delete[] matchcont;
        #elif defined(FourthOrderSpace)
            delete[] match_ouTTercont;
            delete[] match_ouTercont;
            delete[] match_ouTTer;
            delete[] match_ouTer;
            delete[] checkIP_ouTTercont;
            delete[] checkIP_ouTercont;
        #endif
    }

    double** allocate2DArray(int dim1, int dim2) {
        double** array = new double*[dim1];
        for (int i = 0; i < dim1; ++i) {
            array[i] = new double[dim2];
        }
        return array;
    }

    int** allocate2DArrayInt(int dim1, int dim2) {
        int** array = new int*[dim1];
        for (int i = 0; i < dim1; ++i) {
            array[i] = new int[dim2];
        }
        return array;
    }

    double*** allocate3DArray(int dim1, int dim2, int dim3) {
        double*** array = new double**[dim1];
        for (int i = 0; i < dim1; ++i) {
            array[i] = new double*[dim2];
            for (int j = 0; j < dim2; ++j) {
                array[i][j] = new double[dim3];
            }
        }
        return array;
    }
    double dt;
    double dx;
    double dy;
    int nx;
    int ny;
    double TIMEEND;
    double KX;
    double KY;
    double TIMETS;
    double GRAVITY;
    int OUTPUTSTEP;
    double ADDFORCINGTIME;
    std::string OUTPUTPATH;

    double d2t;

    patch csswm[6];

    double **sqrtG;
    double **gamma;
    double ***gLower;
    double ***gUpper;
    double **alpha2D;
    double **beta2D;

    double *sqrtGcont;
    double *gammacont;
    double *gLowercont;
    double *gUppercont;
    double *alpha2Dcont;
    double *beta2Dcont;

    #if defined(EquatorialWave)
        bool status_add_forcing = true;
    #endif

    #if defined(SecondOrderSpace)
        int **match;
        int **checkIP;

        int *matchcont;
        int *checkIPcont;
    #elif defined(FourthOrderSpace)
        int **match_ouTTer;
        int **match_ouTer;
        int **checkIP_ouTTer;
        int **checkIP_ouTer;

        int *match_ouTTercont;
        int *match_ouTercont;
        int *checkIP_ouTTercont;
        int *checkIP_ouTercont;
    #endif


    // ***********************************************************************************
    // In construction.cpp
    static void get_gUpper(double ans[4], double alpha, double beta);
    static void get_gLower(double ans[4], double alpha, double beta);
    static void get_A(double ans[4], int p, double alpha, double beta);
    static void get_IA(double ans[4], int p, double alpha, double beta);
    // ***********************************************************************************

    // ***********************************************************************************
    // In transform.cpp
    static double Cube2Sphere_U(CSSWM &, int, int, int);
    static double Cube2Sphere_V(CSSWM &, int, int, int);
    static double Sphere2Cube_U(CSSWM &, int, int, int);
    static double Sphere2Cube_V(CSSWM &, int, int, int);
    static double Cube2Cube_U(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    static double Cube2Cube_V(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    static double Cube2Cube_BV2AU(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    static double Cube2Cube_BU2AV(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    static double Cube2Cube_U_2(double gLower[4], double IA[4], double A[4], double gUpper[4], double u, double v);
    static double Cube2Cube_V_2(double gLower[4], double IA[4], double A[4], double gUpper[4], double u, double v);
    static void matrixMul(double firstMatrix[4], double secondMatrix[4], double mult[2][2]);
    void Cube2Cube_matrix();
    // ***********************************************************************************

    // ***********************************************************************************
    // In bp_h.cpp
    static void BP_h(CSSWM &);
    #if defined(Mountain)
        void BP_hs(CSSWM &);
    #endif
    static double interpolate(double A1, double A2, double V1, double V2, double B);
    // ***********************************************************************************

    // ***********************************************************************************
    // In bp_wind.cpp
    static void BP_wind_convert(CSSWM &);
    static void BP_wind_interpolation(CSSWM &);
    static void BP_wind_interpolation2(CSSWM &);
    // ***********************************************************************************


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
        static double SteadyGeostrophyH(double, double, double GRAVITY);
        static double SteadyGeostrophyU(double, double);
        static double SteadyGeostrophyV(double);
        static double ConvergenceRateH(double, double);
        static double DeformationalFlowH(double, double);
        static double BarotropicH(double, double GRAVITY);
        static double BarotropicHPrime(double, double);

        static double MountainH(double, double, double GRAVITY);
        static double MountainU(double, double);
        static double MountainV(double);

        static double RossbyHaurwitzH(double, double, double GRAVITY);
        static double RossbyHaurwitzU(double, double);
        static double RossbyHaurwitzV(double, double);

        static double EquatorialWaveH(double, double);

        static double simpson(double, double);
        static double func(double);
    };

    class Outputs {
    public:
        static void grid(CSSWM &);
        static void h(int, CSSWM &);
        static void u(int, CSSWM &);
        static void v(int, CSSWM &);

        static void grid_nc(CSSWM &);
        static void huv_nc(int, CSSWM &);

        static void create_all_directory(CSSWM &);

    private:
        static void create_directory(std::string);
    };

    class Iteration {
    public:
        static void TimeMarching(CSSWM &);
        #if defined(SecondOrderSpace)
            static void ph_pt_2(CSSWM &);
            static void pu_pt_2(CSSWM &);
            static void pv_pt_2(CSSWM &);
        #elif defined(FourthOrderSpace) 
            static void ph_pt_4(CSSWM &);
            static void pu_pt_4(CSSWM &);
            static void pv_pt_4(CSSWM &);
        #endif
        static void nextTimeStep(CSSWM &);
    };

    class NumericalProcess {
    public:
        static void DiffusionAll(CSSWM &);
        static void timeFilterAll(CSSWM &);
    };

private:
    void Construct_gamma_sqrtG_GUpper(double **alpha2D, double **beta2D, double **gamma, double **sqrtG, double ***gUpper, double ***gLower);
    void Construct_p0123_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double **gamma, double **lon, double **lat, double **lon_original, double **x, double **y, double ***A, double ***IA);
    void Construct_p4_lonlat_xy_AIA(double **alpha2D, double **beta2D, double **gamma, double **lon, double **lat, double **lon_original, double **x, double **y, double ***A, double ***IA);
    void Construct_p5_lonlat_xy_AIA(double **alpha2D, double **beta2D, double **gamma, double **lon, double **lat, double **lon_original, double **x, double **y, double ***A, double ***IA);
    #if defined(SecondOrderSpace)
        void initMatch_1point(int **match);
    #elif defined(FourthOrderSpace)
        void initMatch_2point(int **match_ouTer, int **match_ouTTer);
    #endif
};

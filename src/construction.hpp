#include <cmath>
#include <iostream>
#include <iomanip>
#include "define.hpp"


class Config_CSSWM {
public:
    Config_CSSWM(double dt, double dx, double dy)
        : dt(dt), dx(dx), dy(dy) {}
    ~Config_CSSWM() {}

    double dt;
    double dx;
    double dy;
};

class CSSWM {
public:
    CSSWM(const Config_CSSWM &config)
      : dt(config.dt), dx(config.dx), dy(config.dy), 
        #if defined(SecondOrderSpace)
            nx((int) (90/config.dx + 2))
            ny((int) (90/config.dy + 2))
        #elif defined(FourthOrderSpace) 
            nx((int) (90/config.dx + 4)),
            ny((int) (90/config.dy + 4))
        #endif
    {
        allocateMemory();
        initialize();
    }

    ~CSSWM() {
        deallocateMemory();
    }

    void allocateMemory() {
        hp = allocate3DContinuousArray(6, nx, ny, hp_cont);
        h = allocate3DContinuousArray(6, nx, ny, h_cont);
        hm = allocate3DContinuousArray(6, nx, ny, hm_cont);
        up = allocate3DContinuousArray(6, nx, ny, up_cont);
        u = allocate3DContinuousArray(6, nx, ny, u_cont);
        um = allocate3DContinuousArray(6, nx, ny, um_cont);
        vp = allocate3DContinuousArray(6, nx, ny, vp_cont);
        v = allocate3DContinuousArray(6, nx, ny, v_cont);
        vm = allocate3DContinuousArray(6, nx, ny, vm_cont);
        #if defined(Mountain)
            hs = allocate3DContinuousArray(6, nx, ny, hs_cont);
        #endif
        #if defined(EquatorialWave)
            h_forcing = allocate3DContinuousArray(6, nx, ny, h_forcing_cont);
        #endif

        lon = allocate3DContinuousArray(6, nx, ny, lon_cont);
        lat = allocate3DContinuousArray(6, nx, ny, lat_cont);
        lon_original = allocate3DContinuousArray(6, nx, ny, lon_original_cont);
        x = allocate3DContinuousArray(6, nx, ny, x_cont);
        y = allocate3DContinuousArray(6, nx, ny, y_cont);
        A = allocate4DContinuousArray(6, nx, ny, 4, A_cont);
        IA = allocate4DContinuousArray(6, nx, ny, 4, IA_cont);
        #if defined(SecondOrderSpace)
            IP1_L = allocate3DContinuousArray(6, nx, 4, IP1_L_cont);
            IP1_R = allocate3DContinuousArray(6, nx, 4, IP1_R_cont);
            IP1_U = allocate3DContinuousArray(6, nx, 4, IP1_U_cont);
            IP1_D = allocate3DContinuousArray(6, nx, 4, IP1_D_cont);

            match = allocate2DContinuousArray<int>(24, 8, match_cont);
            checkIP = allocate2DContinuousArray<int>(nx, 2, checkIP_cont);
        #elif defined(FourthOrderSpace)
            IP_ouTTer_L = allocate3DContinuousArray(6, nx, 4, IP_ouTTer_L_cont);
            IP_ouTTer_R = allocate3DContinuousArray(6, nx, 4, IP_ouTTer_R_cont);
            IP_ouTTer_U = allocate3DContinuousArray(6, nx, 4, IP_ouTTer_U_cont);
            IP_ouTTer_D = allocate3DContinuousArray(6, nx, 4, IP_ouTTer_D_cont);
            IP_ouTer_L = allocate3DContinuousArray(6, nx, 4, IP_ouTer_L_cont);
            IP_ouTer_R = allocate3DContinuousArray(6, nx, 4, IP_ouTer_R_cont);
            IP_ouTer_U = allocate3DContinuousArray(6, nx, 4, IP_ouTer_U_cont);
            IP_ouTer_D = allocate3DContinuousArray(6, nx, 4, IP_ouTer_D_cont);

            match_ouTTer = allocate2DContinuousArray<int>(24, 8, match_ouTTer_cont);
            match_ouTer = allocate2DContinuousArray<int>(24, 8, match_ouTer_cont);
            checkIP_ouTTer = allocate2DContinuousArray<int>(nx, 2, checkIP_ouTTer_cont);
            checkIP_ouTer = allocate2DContinuousArray<int>(nx, 2, checkIP_ouTer_cont);
        #endif

        sqrtG = allocate2DContinuousArray<double>(nx, ny, sqrtG_cont);
        gamma = allocate2DContinuousArray<double>(nx, ny, gamma_cont);
        gLower = allocate3DContinuousArray(nx, ny, 4, gLower_cont);
        gUpper = allocate3DContinuousArray(nx, ny, 4, gUpper_cont);
        alpha2D = allocate2DContinuousArray<double>(nx, ny, alpha2D_cont);
        beta2D = allocate2DContinuousArray<double>(nx, ny, beta2D_cont);

        return;
    }


    void deallocateMemory() {
        deallocate3DContinuousArray(hp, hp_cont, 6);
        deallocate3DContinuousArray(h, h_cont, 6);
        deallocate3DContinuousArray(hm, hm_cont, 6);
        deallocate3DContinuousArray(up, up_cont, 6);
        deallocate3DContinuousArray(u, u_cont, 6);
        deallocate3DContinuousArray(um, um_cont, 6);
        deallocate3DContinuousArray(vp, vp_cont, 6);
        deallocate3DContinuousArray(v, v_cont, 6);
        deallocate3DContinuousArray(vm, vm_cont, 6);
        #if defined(Mountain)
            deallocate3DContinuousArray(hs, hs_cont, 6);
        #endif
        #if defined(EquatorialWave)
            deallocate3DContinuousArray(h_forcing, h_forcing_cont, 6);
        #endif

        deallocate3DContinuousArray(lon, lon_cont, 6);
        deallocate3DContinuousArray(lat, lat_cont, 6);
        deallocate3DContinuousArray(lon_original, lon_original_cont, 6);
        deallocate3DContinuousArray(x, x_cont, 6);
        deallocate3DContinuousArray(y, y_cont, 6);
        deallocate4DContinuousArray(A, A_cont, 6, nx);
        deallocate4DContinuousArray(IA, IA_cont, 6, nx);

        #if defined(SecondOrderSpace)
            deallocate3DContinuousArray(IP1_L, IP1_L_cont, 6);
            deallocate3DContinuousArray(IP1_R, IP1_R_cont, 6);
            deallocate3DContinuousArray(IP1_U, IP1_U_cont, 6);
            deallocate3DContinuousArray(IP1_D, IP1_D_cont, 6);

            deallocate2DContinuousArray(match, match_cont);
            deallocate2DContinuousArray(checkIP, checkIP_cont);
        #elif defined(FourthOrderSpace)
            deallocate3DContinuousArray(IP_ouTTer_L, IP_ouTTer_L_cont, 6);
            deallocate3DContinuousArray(IP_ouTTer_R, IP_ouTTer_R_cont, 6);
            deallocate3DContinuousArray(IP_ouTTer_U, IP_ouTTer_U_cont, 6);
            deallocate3DContinuousArray(IP_ouTTer_D, IP_ouTTer_D_cont, 6);
            deallocate3DContinuousArray(IP_ouTer_L, IP_ouTer_L_cont, 6);
            deallocate3DContinuousArray(IP_ouTer_R, IP_ouTer_R_cont, 6);
            deallocate3DContinuousArray(IP_ouTer_U, IP_ouTer_U_cont, 6);
            deallocate3DContinuousArray(IP_ouTer_D, IP_ouTer_D_cont, 6);

            deallocate2DContinuousArray(match_ouTTer, match_ouTTer_cont);
            deallocate2DContinuousArray(match_ouTer, match_ouTer_cont);
            deallocate2DContinuousArray(checkIP_ouTTer, checkIP_ouTTer_cont);
            deallocate2DContinuousArray(checkIP_ouTer, checkIP_ouTer_cont);
        #endif

        deallocate2DContinuousArray(sqrtG, sqrtG_cont);
        deallocate2DContinuousArray(gamma, gamma_cont);
        deallocate3DContinuousArray(gLower, gLower_cont, nx);
        deallocate3DContinuousArray(gUpper, gUpper_cont, nx);
        deallocate2DContinuousArray(alpha2D, alpha2D_cont);
        deallocate2DContinuousArray(beta2D, beta2D_cont);
    }

    template <typename T>
    static T** allocate2DContinuousArray(int rows, int cols, T*& contMemory) {
        T** array = new T*[rows];
        contMemory = new T[rows * cols]; // Allocate continuous memory block
        for (int i = 0; i < rows; ++i) {
            array[i] = &contMemory[i * cols]; // Point to segments within continuous block
        }
        return array;
    }

    template <typename T>
    static void deallocate2DContinuousArray(T** array, T* contMemory) {
        if (array != nullptr) {
            delete[] contMemory; // Deallocate the continuous block of memory
            delete[] array;      // Deallocate the array of pointers
        }
    }

    double*** allocate3DContinuousArray(int dim1, int dim2, int dim3, double*& contMemory) {
        double*** array = new double**[dim1];
        contMemory = new double[dim1 * dim2 * dim3]; // Allocate continuous memory block
        for (int i = 0; i < dim1; ++i) {
            array[i] = new double*[dim2];
            for (int j = 0; j < dim2; ++j) {
                array[i][j] = &contMemory[i * dim2 * dim3 + j * dim3]; // Point to segments within continuous block
            }
        }
        return array;
    }

    void deallocate3DContinuousArray(double*** array, double* contMemory, int dim1) {
        if (array != nullptr) {
            delete[] contMemory; // Deallocate the continuous block of memory
            for (int i = 0; i < dim1; ++i) {
                delete[] array[i]; // Deallocate the array of pointers
            }
            delete[] array;      // Deallocate the array of pointers
        }
    }

    double**** allocate4DContinuousArray(int dim1, int dim2, int dim3, int dim4, double*& contMemory) {
        double**** array = new double***[dim1];
        contMemory = new double[dim1 * dim2 * dim3 * dim4]; // Allocate continuous memory block
        for (int i = 0; i < dim1; ++i) {
            array[i] = new double**[dim2];
            for (int j = 0; j < dim2; ++j) {
                array[i][j] = new double*[dim3];
                for (int k = 0; k < dim3; ++k) {
                    array[i][j][k] = &contMemory[i * dim2 * dim3 * dim4 + j * dim3 * dim4 + k * dim4]; // Point to segments within continuous block
                }
            }
        }
        return array;
    }

    void deallocate4DContinuousArray(double**** array, double* contMemory, int dim1, int dim2) {
        if (array != nullptr) {
            delete[] contMemory; // Deallocate the continuous block of memory
            for (int i = 0; i < dim1; ++i) {
                for (int j = 0; j < dim2; ++j) {
                    delete[] array[i][j]; // Deallocate the array of pointers for dim3
                }
                delete[] array[i]; // Deallocate the array of pointers for dim2
            }
            delete[] array; // Deallocate the array of pointers for dim1
        }
    }

    double ***hp;
    double ***h;
    double ***hm;
    double ***up;
    double ***u;
    double ***um;
    double ***vp;
    double ***v;
    double ***vm;
    #if defined(Mountain)
        double ***hs;
    #endif
    #if defined(EquatorialWave)
        double ***h_forcing;
    #endif
    double ***lon;
    double ***lat;
    double ***lon_original;
    double ***x;
    double ***y;
    double ****A;
    double ****IA;
    #if defined(SecondOrderSpace)
        double ***IP1_L;
        double ***IP1_R;
        double ***IP1_U;
        double ***IP1_D;

        int **match;
        int **checkIP;
    #elif defined(FourthOrderSpace)
        double ***IP_ouTTer_L;
        double ***IP_ouTTer_R;
        double ***IP_ouTTer_U;
        double ***IP_ouTTer_D;
        double ***IP_ouTer_L;
        double ***IP_ouTer_R;
        double ***IP_ouTer_U;
        double ***IP_ouTer_D;

        int **match_ouTTer;
        int **match_ouTer;
        int **checkIP_ouTTer;
        int **checkIP_ouTer;
    #endif
    double **sqrtG;
    double **gamma;
    double ***gLower;
    double ***gUpper;
    double **alpha2D;
    double **beta2D;


    double *hp_cont;
    double *h_cont;
    double *hm_cont;
    double *up_cont;
    double *u_cont;
    double *um_cont;
    double *vp_cont;
    double *v_cont;
    double *vm_cont;
    #if defined(Mountain)
        double *hs_cont;
    #endif
    #if defined(EquatorialWave)
        double *h_forcing_cont;
    #endif
    double *lon_cont;
    double *lat_cont;
    double *lon_original_cont;
    double *x_cont;
    double *y_cont;
    double *A_cont;
    double *IA_cont;
    #if defined(SecondOrderSpace)
        double *IP1_L_cont;
        double *IP1_R_cont;
        double *IP1_U_cont;
        double *IP1_D_cont;

        int *match;
        int *checkIP;
    #elif defined(FourthOrderSpace)
        double *IP_ouTTer_L_cont;
        double *IP_ouTTer_R_cont;
        double *IP_ouTTer_U_cont;
        double *IP_ouTTer_D_cont;
        double *IP_ouTer_L_cont;
        double *IP_ouTer_R_cont;
        double *IP_ouTer_U_cont;
        double *IP_ouTer_D_cont;

        int *match_ouTTer_cont;
        int *match_ouTer_cont;
        int *checkIP_ouTTer_cont;
        int *checkIP_ouTer_cont;
    #endif
    double *sqrtG_cont;
    double *gamma_cont;
    double *gLower_cont;
    double *gUpper_cont;
    double *alpha2D_cont;
    double *beta2D_cont;


    double dt;
    double dx;
    double dy;
    int nx;
    int ny;
    #if defined(EquatorialWave)
        bool status_add_forcing = true;
    #endif



    // ***********************************************************************************
    // In construction.cpp
    void get_gUpper(double ans[4], double alpha, double beta);
    void get_gLower(double ans[4], double alpha, double beta);
    void get_A(double ans[4], int p, double alpha, double beta);
    void get_IA(double ans[4], int p, double alpha, double beta);
    void initialize();
    // ***********************************************************************************

    // ***********************************************************************************
    // In transform.cpp
    double Cube2Sphere_U(CSSWM &, int, int, int);
    double Cube2Sphere_V(CSSWM &, int, int, int);
    double Sphere2Cube_U(CSSWM &, int, int, int);
    double Sphere2Cube_V(CSSWM &, int, int, int);
    double Cube2Cube_U(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    double Cube2Cube_V(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    double Cube2Cube_BV2AU(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    double Cube2Cube_BU2AV(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    double Cube2Cube_U_2(double gLower[4], double IA[4], double A[4], double gUpper[4], double u, double v);
    double Cube2Cube_V_2(double gLower[4], double IA[4], double A[4], double gUpper[4], double u, double v);
    void matrixMul(double firstMatrix[4], double secondMatrix[4], double mult[2][2]);
    void Cube2Cube_matrix();
    // ***********************************************************************************

    // ***********************************************************************************
    // In bp_h.cpp
    void BP_h(CSSWM &);
    #if defined(Mountain)
        void BP_hs(CSSWM &);
    #endif
    double interpolate(double A1, double A2, double V1, double V2, double B);
    // ***********************************************************************************

    // ***********************************************************************************
    // In bp_wind.cpp
    void BP_wind_convert(CSSWM &);
    void BP_wind_interpolation(CSSWM &);
    void BP_wind_interpolation2(CSSWM &);
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

        static void create_all_directory();

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
    void Construct_p4_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double **gamma, double **lon, double **lat, double **lon_original, double **x, double **y, double ***A, double ***IA);
    void Construct_p5_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double **gamma, double **lon, double **lat, double **lon_original, double **x, double **y, double ***A, double ***IA);
};

#include "construction.hpp"
#include <netcdf>
#include <fstream>

using std::fstream;
using std::ios;
using std::string;
using namespace netCDF;

void CSSWM::Outputs::create_directory(string directory_name) {
    string str = "mkdir -p " + directory_name;
    const char *command = str.c_str();
    const int dir_err = system(command);
    if (-1 == dir_err) {
        std::cout << "Error on creating directory!\n" << std::endl;
        return;
    }
    return;
}

void CSSWM::Outputs::grid(CSSWM &model) {
    fstream fout[4];
    string dir = model.outputpath + (string) "grids/";
    string grid[4] = {"lon.txt", "lat.txt", "x.txt", "y.txt"};

    for (int i = 0; i < 4; i++) {
        fout[i].open(dir + grid[i], ios::out);
    }

    for (int p = 0; p < 6; p++) {
        for (int j = 1; j < model.ny-1; j++) {
            for (int i = 1; i < model.nx-1; i++) {
                fout[0] << model.lon[p][i][j] << " ";
                fout[1] << model.lat[p][i][j] << " ";
        
                fout[2] << model.x[p][i][j] << " ";
                fout[3] << model.y[p][i][j] << " ";
            }
        }
    }
}

void CSSWM::Outputs::h(int n, CSSWM &model) {
    fstream fouth;
    string hname = model.outputpath + (string) "h/h_" + std::to_string(n) + ".txt";
    fouth.open(hname, std::ios::out);
    for (int p = 0; p < 6; p++) {
        for (int j = 1; j < model.ny-1; j++) {
            for (int i = 1; i < model.nx-1; i++) {
                fouth << model.h[p][i][j] << " ";
            }
        }
    }
    return;
}

void CSSWM::Outputs::u(int n, CSSWM &model) {
    fstream foutu;
    string uname = model.outputpath + (string) "u/u_" + std::to_string(n) + ".txt";
    foutu.open(uname, std::ios::out);

    fstream foutu_lon_lat;
    string u_lon_latname = "../outputs/u_lon_lat/u_lon_lat_" + std::to_string(n) + ".txt";
    foutu_lon_lat.open(u_lon_latname, std::ios::out);
    for (int p = 0; p < 6; p++) {
        for (int j = 1; j < model.ny-1; j++) {
            for (int i = 1; i < model.nx-1; i++) {
                foutu << model.u[p][i][j] << " ";
                foutu_lon_lat << model.Cube2Sphere_U(model, p, i, j) << " ";
            }
        }
    }
    return;
}

void CSSWM::Outputs::v(int n, CSSWM &model) {
    fstream foutv;
    string vname = model.outputpath + (string) "v/v_" + std::to_string(n) + ".txt";
    foutv.open(vname, std::ios::out);

    fstream foutv_lon_lat;
    string v_lon_latname = "../outputs/v_lon_lat/v_lon_lat_" + std::to_string(n) + ".txt";
    foutv_lon_lat.open(v_lon_latname, std::ios::out);
    for (int p = 0; p < 6; p++) {
        for (int j = 1; j < model.ny-1; j++) {
            for (int i = 1; i < model.nx-1; i++) {
                foutv << model.v[p][i][j] << " ";
                foutv_lon_lat << model.Cube2Sphere_V(model, p, i, j) << " ";
            }
        }
    }
    return;
}


void checkErr(int status, int line) {
    if (status != NC_NOERR) {
        std::cerr << "NetCDF error at line " << line << ": " << nc_strerror(status) << std::endl;
        exit(EXIT_FAILURE);
    }
}

void CSSWM::Outputs::grid_nc(CSSWM &model) {
    string ncName = model.outputpath + (string) "nc/" + "grid.nc";

    int ncid, p_dimid, x_dimid, y_dimid;
    int retval;

    int xid, yid, lonid, latid;

    if ((retval = nc_create(ncName.c_str(), NC_CLOBBER, &ncid))) checkErr(retval, __LINE__);

    if ((retval = nc_def_dim(ncid, "p", 6, &p_dimid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_dim(ncid, "x", model.nx, &x_dimid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_dim(ncid, "y", model.ny, &y_dimid))) checkErr(retval, __LINE__);

    int dimids[3] = {p_dimid, x_dimid, y_dimid};
    // int dimx1d[1] = {x_dimid};
    // int dimz1d[1] = {z_dimid};
    if ((retval = nc_def_var(ncid, "x_local", NC_DOUBLE, 3, dimids, &xid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "y_local", NC_DOUBLE, 3, dimids, &yid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "lon_sphere", NC_DOUBLE, 3, dimids, &lonid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "lat_sphere", NC_DOUBLE, 3, dimids, &latid))) checkErr(retval, __LINE__);

    if ((retval = nc_enddef(ncid))) checkErr(retval, __LINE__);

    if ((retval = nc_put_var_double(ncid, xid, model.x_cont))) checkErr(retval, __LINE__);
    if ((retval = nc_put_var_double(ncid, yid, model.y_cont))) checkErr(retval, __LINE__);
    if ((retval = nc_put_var_double(ncid, lonid, model.lon_cont))) checkErr(retval, __LINE__);
    if ((retval = nc_put_var_double(ncid, latid, model.lat_cont))) checkErr(retval, __LINE__);


    if ((retval = nc_close(ncid))) checkErr(retval, __LINE__);
}


void CSSWM::Outputs::huv_nc(int n, CSSWM &model) {
    string ncName = model.outputpath + (string) "nc/" + std::to_string(n) + ".nc";

    int ncid, p_dimid, x_dimid, y_dimid;
    int retval;

    int hid, uid, vid;
    // int ulonlatid, vlonlatid;

    if ((retval = nc_create(ncName.c_str(), NC_CLOBBER, &ncid))) checkErr(retval, __LINE__);

    if ((retval = nc_def_dim(ncid, "p", 6, &p_dimid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_dim(ncid, "x", model.nx, &x_dimid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_dim(ncid, "y", model.ny, &y_dimid))) checkErr(retval, __LINE__);

    int dimids[3] = {p_dimid, x_dimid, y_dimid};
    if ((retval = nc_def_var(ncid, "h", NC_DOUBLE, 3, dimids, &hid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "u", NC_DOUBLE, 3, dimids, &uid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "v", NC_DOUBLE, 3, dimids, &vid))) checkErr(retval, __LINE__);
    // if ((retval = nc_def_var(ncid, "u_lonlat", NC_DOUBLE, 3, dimids, &ulonlatid))) checkErr(retval, __LINE__);
    // if ((retval = nc_def_var(ncid, "v_lonlat", NC_DOUBLE, 3, dimids, &vlonlatid))) checkErr(retval, __LINE__);

    if ((retval = nc_enddef(ncid))) checkErr(retval, __LINE__);

    if ((retval = nc_put_var_double(ncid, hid, model.h_cont))) checkErr(retval, __LINE__);
    if ((retval = nc_put_var_double(ncid, uid, model.u_cont))) checkErr(retval, __LINE__);
    if ((retval = nc_put_var_double(ncid, vid, model.v_cont))) checkErr(retval, __LINE__);
    // if ((retval = nc_put_var_double(ncid, ulonlatid, model.u_cont))) checkErr(retval, __LINE__);
    // if ((retval = nc_put_var_double(ncid, vlonlatid, model.v_cont))) checkErr(retval, __LINE__);


    if ((retval = nc_close(ncid))) checkErr(retval, __LINE__);
}

void CSSWM::Outputs::create_all_directory(std::string outputpath) {
    // data directory
    #ifdef TXTOUTPUT
        create_directory(model.outputpath + (string) "grids");
        create_directory(model.outputpath + (string) "h");
        create_directory(model.outputpath + (string) "u");
        create_directory(model.outputpath + (string) "u_lon_lat");
        create_directory(model.outputpath + (string) "v");
        create_directory(model.outputpath + (string) "v_lon_lat");
    #endif
    #ifdef NCOUTPUT
        create_directory(outputpath + (string) "nc");
    #endif
}
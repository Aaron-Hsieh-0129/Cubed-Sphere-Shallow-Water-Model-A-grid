#include "construction.hpp"
#include <vector>
#include <netcdf>
#include <fstream>

using std::fstream;
using std::ios;
using std::string;
using std::vector;
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
    string dir = OUTPUTPATH + (string) "grids/";
    string grid[4] = {"lon.txt", "lat.txt", "x.txt", "y.txt"};

    for (int i = 0; i < 4; i++) {
        fout[i].open(dir + grid[i], ios::out);
    }

    for (int p = 0; p < 6; p++) {
        for (int j = 1; j < NY-1; j++) {
            for (int i = 1; i < NX-1; i++) {
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
    string hname = OUTPUTPATH + (string) "h/h_" + std::to_string(n) + ".txt";
    fouth.open(hname, std::ios::out);
    for (int p = 0; p < 6; p++) {
        for (int j = 1; j < NY-1; j++) {
            for (int i = 1; i < NX-1; i++) {
                fouth << model.h[p][i][j] << " ";
            }
        }
    }
    return;
}

void CSSWM::Outputs::u(int n, CSSWM &model) {
    fstream foutu;
    string uname = OUTPUTPATH + (string) "u/u_" + std::to_string(n) + ".txt";
    foutu.open(uname, std::ios::out);

    fstream foutu_lon_lat;
    string u_lon_latname = "../outputs/u_lon_lat/u_lon_lat_" + std::to_string(n) + ".txt";
    foutu_lon_lat.open(u_lon_latname, std::ios::out);
    for (int p = 0; p < 6; p++) {
        for (int j = 1; j < NY-1; j++) {
            for (int i = 1; i < NX-1; i++) {
                foutu << model.u[p][i][j] << " ";
                foutu_lon_lat << model.Cube2Sphere_U(model, p, i, j) << " ";
            }
        }
    }
    return;
}

void CSSWM::Outputs::v(int n, CSSWM &model) {
    fstream foutv;
    string vname = OUTPUTPATH + (string) "v/v_" + std::to_string(n) + ".txt";
    foutv.open(vname, std::ios::out);

    fstream foutv_lon_lat;
    string v_lon_latname = "../outputs/v_lon_lat/v_lon_lat_" + std::to_string(n) + ".txt";
    foutv_lon_lat.open(v_lon_latname, std::ios::out);
    for (int p = 0; p < 6; p++) {
        for (int j = 1; j < NY-1; j++) {
            for (int i = 1; i < NX-1; i++) {
                foutv << model.v[p][i][j] << " ";
                foutv_lon_lat << model.Cube2Sphere_V(model, p, i, j) << " ";
            }
        }
    }
    return;
}

// void CSSWM::Outputs::grid_nc(CSSWM &model) {
//     string dir = OUTPUTPATH + (string) "nc/";

//     NcFile dataFile(dir + "grid.nc", NcFile::replace);       
//     // Create netCDF dimensions
//     NcDim p = dataFile.addDim("p", 6);
//     NcDim xDim = dataFile.addDim("x", NX);
//     NcDim yDim = dataFile.addDim("y", NY);
//     NcDim lonDim = dataFile.addDim("lon", NX);
//     NcDim latDim = dataFile.addDim("lat", NY);

//     vector<NcDim> xyDim, lonlatDim;
//     xyDim.push_back(p);
//     xyDim.push_back(xDim);
//     xyDim.push_back(yDim);

//     lonlatDim.push_back(p);
//     lonlatDim.push_back(lonDim);
//     lonlatDim.push_back(latDim);

//     NcVar x = dataFile.addVar("x_local", ncDouble, xyDim);
//     NcVar y = dataFile.addVar("y_local", ncDouble, xyDim);
//     NcVar lon = dataFile.addVar("lon_sphere", ncDouble, lonlatDim);
//     NcVar lat = dataFile.addVar("lat_sphere", ncDouble, lonlatDim);
//     NcVar A = dataFile.addVar("area_sphere_coeff", ncDouble, lonlatDim);
//     #if defined(Mountain)
//         NcVar hs = dataFile.addVar("hs", ncDouble, xyDim);
//     #endif

//     double area[6][NX][NY];
//     for (int p = 0; p < 6; p++) {
//         for (int i = 0; i < NX; i++) {
//             for (int j = 0; j < NY; j++) {
//                 area[p][i][j] = model.sqrtG[i][j];
//             }
//         }
//     }
    
//     vector<size_t> startp, countp;
//     startp.push_back(0);
//     startp.push_back(0);
//     startp.push_back(0);
//     countp.push_back(1);
//     countp.push_back(NX);
//     countp.push_back(NY);

//     for (int p = 0; p < 6; p++) {
//         startp[0] = p;
//         x.putVar(startp, countp, model.x[p]);
//         y.putVar(startp, countp, model.y[p]);
//         lon.putVar(startp, countp, model.lon[p]);
//         lat.putVar(startp, countp, model.lat[p]);
//         A.putVar(startp, countp, area);
//         #if defined(Mountain)
//             hs.putVar(startp, countp, model.csswm[p].hs);
//         #endif
//     }
// }
void checkErr(int status, int line) {
    if (status != NC_NOERR) {
        std::cerr << "NetCDF error at line " << line << ": " << nc_strerror(status) << std::endl;
        exit(EXIT_FAILURE);
    }
}

void CSSWM::Outputs::grid_nc(CSSWM &model) {
    string ncName = OUTPUTPATH + (string) "nc/" + "grid.nc";

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

// void CSSWM::Outputs::huv_nc(int n, CSSWM &model) {
//     string dir = OUTPUTPATH + (string) "nc/";

//     NcFile dataFile(dir + std::to_string(n) + ".nc", NcFile::replace);       
//     // Create netCDF dimensions
//     NcDim p = dataFile.addDim("p", 6);
//     NcDim xDim = dataFile.addDim("x", NX);
//     NcDim yDim = dataFile.addDim("y", NY);
//     NcDim lonDim = dataFile.addDim("lon", NX);
//     NcDim latDim = dataFile.addDim("lat", NY);

//     vector<NcDim> xyDim, lonlatDim;
//     xyDim.push_back(p);
//     xyDim.push_back(xDim);
//     xyDim.push_back(yDim);

//     lonlatDim.push_back(p);
//     lonlatDim.push_back(lonDim);
//     lonlatDim.push_back(latDim);

//     NcVar h = dataFile.addVar("h", ncDouble, xyDim);
//     NcVar u = dataFile.addVar("u", ncDouble, xyDim);
//     NcVar v = dataFile.addVar("v", ncDouble, xyDim);

//     NcVar ulonlat = dataFile.addVar("u_lonlat", ncDouble, lonlatDim);
//     NcVar vlonlat = dataFile.addVar("v_lonlat", ncDouble, lonlatDim);
//     double u_lon_lat[6][NX][NY], v_lon_lat[6][NX][NY];
//     for (int p = 0; p < 6; p++) {
//         for (int j = 0; j < NY; j++) {
//             for (int i = 0; i < NX; i++) {
//                 u_lon_lat[p][i][j] = model.Cube2Sphere_U(model, p, i, j);
//                 v_lon_lat[p][i][j] = model.Cube2Sphere_V(model, p, i, j);
//             }
//         }
//     }

//     vector<size_t> startp, countp;
//     startp.push_back(0);
//     startp.push_back(0);
//     startp.push_back(0);
//     countp.push_back(1);
//     countp.push_back(NX);
//     countp.push_back(NY);

//     for (int p = 0; p < 6; p++) {
//         startp[0] = p;
//         h.putVar(startp, countp, model.h[p]);
//         u.putVar(startp, countp, model.u[p]);
//         v.putVar(startp, countp, model.v[p]);

//         ulonlat.putVar(startp, countp, u_lon_lat[p]);
//         vlonlat.putVar(startp, countp, v_lon_lat[p]);
//     }
// }


void CSSWM::Outputs::huv_nc(int n, CSSWM &model) {
    string ncName = OUTPUTPATH + (string) "nc/" + std::to_string(n) + ".nc";

    int ncid, p_dimid, x_dimid, y_dimid;
    int retval;

    int hid, uid, vid, ulonlatid, vlonlatid;

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

void CSSWM::Outputs::create_all_directory() {
    // data directory
    #ifdef TXTOUTPUT
        create_directory(OUTPUTPATH + (string) "grids");
        create_directory(OUTPUTPATH + (string) "h");
        create_directory(OUTPUTPATH + (string) "u");
        create_directory(OUTPUTPATH + (string) "u_lon_lat");
        create_directory(OUTPUTPATH + (string) "v");
        create_directory(OUTPUTPATH + (string) "v_lon_lat");
    #endif
    #ifdef NCOUTPUT
        create_directory(OUTPUTPATH + (string) "nc");
    #endif
}
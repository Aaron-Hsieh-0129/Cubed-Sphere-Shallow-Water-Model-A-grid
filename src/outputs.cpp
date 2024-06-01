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
                fout[0] << model.csswm[p].lon[i][j] << " ";
                fout[1] << model.csswm[p].lat[i][j] << " ";
        
                fout[2] << model.csswm[p].x[i][j] << " ";
                fout[3] << model.csswm[p].y[i][j] << " ";
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
                fouth << model.csswm[p].h[i][j] << " ";
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
                foutu << model.csswm[p].u[i][j] << " ";
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
                foutv << model.csswm[p].v[i][j] << " ";
                foutv_lon_lat << model.Cube2Sphere_V(model, p, i, j) << " ";
            }
        }
    }
    return;
}

void CSSWM::Outputs::grid_nc(CSSWM &model) {
    string dir = OUTPUTPATH + (string) "nc/";

    NcFile dataFile(dir + "grid.nc", NcFile::replace);       
    // Create netCDF dimensions
    NcDim p = dataFile.addDim("p", 6);
    NcDim xDim = dataFile.addDim("x", NX);
    NcDim yDim = dataFile.addDim("y", NY);
    NcDim lonDim = dataFile.addDim("lon", NX);
    NcDim latDim = dataFile.addDim("lat", NY);

    vector<NcDim> xyDim, lonlatDim;
    xyDim.push_back(p);
    xyDim.push_back(xDim);
    xyDim.push_back(yDim);

    lonlatDim.push_back(p);
    lonlatDim.push_back(lonDim);
    lonlatDim.push_back(latDim);

    NcVar x = dataFile.addVar("x_local", ncDouble, xyDim);
    NcVar y = dataFile.addVar("y_local", ncDouble, xyDim);
    NcVar lon = dataFile.addVar("lon_sphere", ncDouble, lonlatDim);
    NcVar lat = dataFile.addVar("lat_sphere", ncDouble, lonlatDim);
    NcVar A = dataFile.addVar("area_sphere_coeff", ncDouble, lonlatDim);
    #if defined(Mountain)
        NcVar hs = dataFile.addVar("hs", ncDouble, xyDim);
    #endif

    double area[6][NX][NY];
    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                area[p][i][j] = model.sqrtG[i][j];
            }
        }
    }
    
    vector<size_t> startp, countp;
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);
    countp.push_back(1);
    countp.push_back(NX);
    countp.push_back(NY);

    for (int p = 0; p < 6; p++) {
        startp[0] = p;
        x.putVar(startp, countp, model.csswm[p].x);
        y.putVar(startp, countp, model.csswm[p].y);
        lon.putVar(startp, countp, model.csswm[p].lon);
        lat.putVar(startp, countp, model.csswm[p].lat);
        A.putVar(startp, countp, area);
        #if defined(Mountain)
            hs.putVar(startp, countp, model.csswm[p].hs);
        #endif
    }
}

void CSSWM::Outputs::huv_nc(int n, CSSWM &model) {
    string dir = OUTPUTPATH + (string) "nc/";

    NcFile dataFile(dir + std::to_string(n) + ".nc", NcFile::replace);       
    // Create netCDF dimensions
    NcDim p = dataFile.addDim("p", 6);
    NcDim xDim = dataFile.addDim("x", NX);
    NcDim yDim = dataFile.addDim("y", NY);
    NcDim lonDim = dataFile.addDim("lon", NX);
    NcDim latDim = dataFile.addDim("lat", NY);

    vector<NcDim> xyDim, lonlatDim;
    xyDim.push_back(p);
    xyDim.push_back(xDim);
    xyDim.push_back(yDim);

    lonlatDim.push_back(p);
    lonlatDim.push_back(lonDim);
    lonlatDim.push_back(latDim);

    NcVar h = dataFile.addVar("h", ncDouble, xyDim);
    NcVar u = dataFile.addVar("u", ncDouble, xyDim);
    NcVar v = dataFile.addVar("v", ncDouble, xyDim);

    NcVar ulonlat = dataFile.addVar("u_lonlat", ncDouble, lonlatDim);
    NcVar vlonlat = dataFile.addVar("v_lonlat", ncDouble, lonlatDim);
    double u_lon_lat[6][NX][NY], v_lon_lat[6][NX][NY];
    for (int p = 0; p < 6; p++) {
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                u_lon_lat[p][i][j] = model.Cube2Sphere_U(model, p, i, j);
                v_lon_lat[p][i][j] = model.Cube2Sphere_V(model, p, i, j);
            }
        }
    }

    vector<size_t> startp, countp;
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);
    countp.push_back(1);
    countp.push_back(NX);
    countp.push_back(NY);

    for (int p = 0; p < 6; p++) {
        startp[0] = p;
        h.putVar(startp, countp, model.csswm[p].h);
        u.putVar(startp, countp, model.csswm[p].u);
        v.putVar(startp, countp, model.csswm[p].v);

        ulonlat.putVar(startp, countp, u_lon_lat[p]);
        vlonlat.putVar(startp, countp, v_lon_lat[p]);
    }
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
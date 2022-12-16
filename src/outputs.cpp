#include "outputs.hpp"

using std::fstream;
using std::ios;
using std::string;

void Outputs::create_directory(string directory_name) {
    string str = "mkdir -p " + directory_name;
    const char *command = str.c_str();
    const int dir_err = system(command);
    if (-1 == dir_err) {
        std::cout << "Error on creating directory!\n" << std::endl;
        return;
    }
    return;
}

void Outputs::output_parameter(CSSWM &model) {
    create_directory("../outputs/grids");
    create_directory("../graphs/grids");

    fstream fout[4];
    string dir = "../outputs/grids/";
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

void Outputs::output_h(int n, CSSWM &model) {
    create_directory("../outputs/h");
    create_directory("../graphs/h/curvilinear");
    create_directory("../graphs/h/sphere");

    fstream fouth;
    string hname = "../outputs/h/h_" + std::to_string(n) + ".txt";
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

void Outputs::output_u(int n, CSSWM &model) {
    create_directory("../outputs/u");
    create_directory("../outputs/u_lon_lat");
    
    fstream foutu;
    string uname = "../outputs/u/u_" + std::to_string(n) + ".txt";
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

void Outputs::output_v(int n, CSSWM &model) {
    create_directory("../outputs/v");
    create_directory("../outputs/v_lon_lat");

    fstream foutv;
    string vname = "../outputs/v/v_" + std::to_string(n) + ".txt";
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
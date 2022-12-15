#include "constrcution.hpp"


void CSSWM::matrixMul(double firstMatrix[4], double secondMatrix[4], double mult[2][2]) {
	double A[2][2], B[2][2];

    // Init
    int count = 0;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			mult[i][j] = 0;
            A[i][j] = firstMatrix[count];
            B[i][j] = secondMatrix[count];
            count++;
		}
	}
    // multiplication
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				mult[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

double CSSWM::Cube2Sphere_U(CSSWM &model, int p, int i, int j) {
    double mult[2][2];
    model.matrixMul(model.gUpper[i][j], model.csswm[p].A[i][j], mult);
    return mult[0][0] * model.csswm[p].u[i][j] + mult[0][1] * model.csswm[p].v[i][j];
}

double CSSWM::Cube2Sphere_V(CSSWM &model, int p, int i, int j) {
    double mult[2][2];
    matrixMul(model.csswm[p].A[i][j], model.gUpper[i][j], mult);
    return mult[1][0] * model.csswm[p].u[i][j] + mult[1][1] * model.csswm[p].v[i][j];
}

double CSSWM::Sphere2Cube_U(CSSWM &model, int p, int i, int j) {
    double mult[2][2];
    matrixMul(model.gLower[i][j], model.csswm[p].IA[i][j], mult);
    return mult[0][0] * model.csswm[p].u[i][j] + mult[0][1] * model.csswm[p].v[i][j];
}

double CSSWM::Sphere2Cube_V(CSSWM &model, int p, int i, int j) {
    double mult[2][2];
    matrixMul(model.gLower[i][j], model.csswm[p].IA[i][j], mult);
    return mult[1][0] * model.csswm[p].u[i][j] + mult[1][1] * model.csswm[p].v[i][j];
}

double CSSWM::Cube2Cube_U(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2) {
    // p1 is other patch and p2 is the patch who needs other's information
    double mult[2][2], A[2][2], B[2][2];
    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            mult[i][j] = A[i][j] = B[i][j] = 0.;
        }
    }

    matrixMul(model.gLower[i1][j1], model.csswm[p1].IA[i1][j1], A);

    matrixMul(model.csswm[p2].A[i2][j2], model.gUpper[i2][j2], B);
    

    // multiply A & B
    for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				mult[i][j] += A[i][k] * B[k][j];
			}
		}
	}
    
    return mult[0][0] * model.csswm[p2].up[i2][j2] + mult[0][1] * model.csswm[p2].vp[i2][j2];
}

double CSSWM::Cube2Cube_V(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2) {
    // p1 is other patch and p2 is the patch who needs other's information
    double mult[2][2], A[2][2], B[2][2];
    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            mult[i][j] = A[i][j] = B[i][j] = 0.;
        }
    }

    matrixMul(model.gLower[i1][j1], model.csswm[p1].IA[i1][j1], A);

    matrixMul(model.csswm[p2].A[i2][j2], model.gUpper[i2][j2], B);
    

    // multiply A & B
    for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				mult[i][j] += A[i][k] * B[k][j];
			}
		}
	}

    return mult[1][0] * model.csswm[p2].up[i2][j2] + mult[1][1] * model.csswm[p2].vp[i2][j2];
}

double CSSWM::Cube2Cube_BV2AU(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2) {
    // p1 is other patch and p2 is the patch who needs other's information
    double mult[2][2], A[2][2], B[2][2];
    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            mult[i][j] = A[i][j] = B[i][j] = 0.;
        }
    }

    matrixMul(model.gLower[i1][j1], model.csswm[p1].IA[i1][j1], A);

    matrixMul(model.csswm[p2].A[i2][j2], model.gUpper[i2][j2], B);
    

    // multiply A & B
    for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				mult[i][j] += A[i][k] * B[k][j];
			}
		}
	}

    return mult[0][0] * model.csswm[p2].up[i2][j2] + mult[0][1] * model.csswm[p2].vp[i2][j2];
}

double CSSWM::Cube2Cube_BU2AV(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2) {
    // p1 is other patch and p2 is the patch who needs other's information
    double mult[2][2], A[2][2], B[2][2];
    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            mult[i][j] = A[i][j] = B[i][j] = 0.;
        }
    }

    matrixMul(model.gLower[i1][j1], model.csswm[p1].IA[i1][j1], A);

    matrixMul(model.csswm[p2].A[i2][j2], model.gUpper[i2][j2], B);
    

    // multiply A & B
    for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				mult[i][j] += A[i][k] * B[k][j];
			}
		}
	}

    return mult[1][0] * model.csswm[p2].up[i2][j2] + mult[1][1] * model.csswm[p2].vp[i2][j2];
}
#include "construction.hpp"


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

    // test
    // model.csswm[p2].up[i2][j2] = -10.;
    // model.csswm[p2].vp[i2][j2] = 0.;
    // std::cout << "u: " << mult[0][0] * model.csswm[p2].up[i2][j2] + mult[0][1] * model.csswm[p2].vp[i2][j2] << std::endl;
    // std::cout << "v: " << mult[1][0] * model.csswm[p2].up[i2][j2] + mult[1][1] * model.csswm[p2].vp[i2][j2] << std::endl;

    return mult[1][0] * model.csswm[p2].up[i2][j2] + mult[1][1] * model.csswm[p2].vp[i2][j2];
}

double CSSWM::Cube2Cube_U_2(double gLower[4], double IA[4], double A[4], double gUpper[4], double u, double v) {
    double mult[2][2], tmp1[2][2], tmp2[2][2];

    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            mult[i][j] = tmp1[i][j] = tmp2[i][j] = 0.;
        }
    }

    matrixMul(gLower, IA, tmp1);

    matrixMul(A, gUpper, tmp2);

    for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				mult[i][j] += tmp1[i][k] * tmp2[k][j];
			}
		}
	}

    return mult[0][0] * u + mult[0][1] * v;
}

double CSSWM::Cube2Cube_V_2(double gLower[4], double IA[4], double A[4], double gUpper[4], double u, double v) {
    double mult[2][2], tmp1[2][2], tmp2[2][2];

    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            mult[i][j] = tmp1[i][j] = tmp2[i][j] = 0.;
        }
    }

    matrixMul(gLower, IA, tmp1);

    matrixMul(A, gUpper, tmp2);

    for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				mult[i][j] += tmp1[i][k] * tmp2[k][j];
			}
		}
	}

    return mult[1][0] * u + mult[1][1] * v;
}


void CSSWM::Cube2Cube_matrix() {
    // Construct the patch to patch transformation matrix 
    int p1, p2, i1, j1, i2, j2, reversed, lonlat;
    double alpha_tmp, beta_tmp, alpha_B_tmp, beta_B_tmp;
    double tmp1[2][2], tmp2[2][2], mult[2][2];
    int count = 0;
    double A1, A2, B;
    double gLower_tmp[4], gUpper_tmp[4], A_tmp[4], IA_tmp[4];
    int I1 = -999, I2_1 = -999, I2_2 = -999, J1 = -999, J2_1 = -999, J2_2 = -999;
    for (int pp = 0; pp < 24; pp++) {
        #if defined(SecondOrderSpace)
            p1 = match[pp][0], p2 = match[pp][1], i1 = match[pp][2], j1 = match[pp][3], i2 = match[pp][4], j2 = match[pp][5], reversed = match[pp][6], lonlat = match[pp][7];
        #elif defined(FourthOrderSpace)
            p1 = match_ouTTer[pp][0], p2 = match_ouTTer[pp][1], i1 = match_ouTTer[pp][2], j1 = match_ouTTer[pp][3], i2 = match_ouTTer[pp][4], j2 = match_ouTTer[pp][5], reversed = match_ouTTer[pp][6], lonlat = match_ouTTer[pp][7];
        #endif
        for (int idx = 0; idx < NX; idx++) {
            I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
            #if defined(SecondOrderSpace)
                I2_1 = i2 == -1 ? reversed ? checkIP[NX-1-idx][0] : checkIP[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? checkIP[NY-1-idx][0] : checkIP[idx][0] : j2;
                I2_2 = i2 == -1 ? reversed ? checkIP[NX-1-idx][1] : checkIP[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? checkIP[NY-1-idx][1] : checkIP[idx][1] : j2;
            #elif defined(FourthOrderSpace)
                I2_1 = i2 == -1 ? reversed ? checkIP_ouTTer[NX-1-idx][0] : checkIP_ouTTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? checkIP_ouTTer[NY-1-idx][0] : checkIP_ouTTer[idx][0] : j2;
                I2_2 = i2 == -1 ? reversed ? checkIP_ouTTer[NX-1-idx][1] : checkIP_ouTTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? checkIP_ouTTer[NY-1-idx][1] : checkIP_ouTTer[idx][1] : j2;
            #endif

            if (lonlat == 0) {
                B = csswm[p1].lat[I1][J1];
                A1 = csswm[p2].lat[I2_1][J2_1], A2 = csswm[p2].lat[I2_2][J2_2];
            }
            else {
                B = csswm[p1].lon[I1][J1];
                A1 = csswm[p2].lon[I2_1][J2_1], A2 = csswm[p2].lon[I2_2][J2_2];

                if (A1 > A2 && (p1 == 0 || p2 == 0))  A2 += 2 * M_PI;
                if (A1 > B && B < A2) B += 2 * M_PI;
            }

            alpha_tmp = alpha2D[I1][J1];
            beta_tmp = beta2D[I1][J1];

            alpha_B_tmp = interpolate(A1, A2, alpha2D[I2_1][J2_1], alpha2D[I2_2][J2_2], B);
            beta_B_tmp = interpolate(A1, A2, beta2D[I2_1][J2_1], beta2D[I2_2][J2_2], B);

            get_gLower(gLower_tmp, alpha_tmp, beta_tmp);
            get_IA(IA_tmp, p1, alpha_tmp, beta_tmp);
            get_A(A_tmp, p2, alpha_B_tmp, beta_B_tmp);
            get_gUpper(gUpper_tmp, alpha_B_tmp, beta_B_tmp);

            matrixMul(gLower_tmp, IA_tmp, tmp1);
            matrixMul(A_tmp, gUpper_tmp, tmp2);
            double tmp1_4[4], tmp2_4[4];
            count = 0;
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    tmp1_4[count] = tmp1[i][j];
                    tmp2_4[count] = tmp2[i][j];
                    count++;
                }
            }
            matrixMul(tmp1_4, tmp2_4, mult);
            
            if (i1 == -1) {
                if (j1 == 0) {
                    count = 0;
                    for (int i = 0; i < 2; i++) {
                        for (int j = 0; j < 2; j++) {
                            #if defined(SecondOrderSpace)
                                csswm[p1].IP1_D[I1][count] = mult[i][j];
                            #elif defined(FourthOrderSpace)
                                csswm[p1].IP_ouTTer_D[I1][count] = mult[i][j];
                            #endif
                            count++;
                        }
                    }
                }
                else if (j1 == NY-1) {
                    count = 0;
                    for (int i = 0; i < 2; i++) {
                        for (int j = 0; j < 2; j++) {
                            #if defined(SecondOrderSpace)
                                csswm[p1].IP1_U[I1][count] = mult[i][j];
                            #elif defined(FourthOrderSpace)
                                csswm[p1].IP_ouTTer_U[I1][count] = mult[i][j];
                            #endif
                            count++;
                        }
                    }
                }
            }
            else if (j1 == -1) {
                if (i1 == 0) {
                    count = 0;
                    for (int i = 0; i < 2; i++) {
                        for (int j = 0; j < 2; j++) {
                            #if defined(SecondOrderSpace)
                                csswm[p1].IP1_L[J1][count] = mult[i][j];
                            #elif defined(FourthOrderSpace)
                                csswm[p1].IP_ouTTer_L[J1][count] = mult[i][j];
                            #endif
                            count++;
                        }
                    }
                }
                else if (i1 == NY-1) {
                    count = 0;
                    for (int i = 0; i < 2; i++) {
                        for (int j = 0; j < 2; j++) {
                            #if defined(SecondOrderSpace)
                                csswm[p1].IP1_R[J1][count] = mult[i][j];
                            #elif defined(FourthOrderSpace)
                                csswm[p1].IP_ouTTer_R[J1][count] = mult[i][j];
                            #endif
                            count++;
                        }
                    }
                }
            }
        }
    }


    #if defined(FourthOrderSpace)
        // Construct the patch to patch transformation matrix 
        for (int pp = 0; pp < 24; pp++) {
            p1 = match_ouTer[pp][0], p2 = match_ouTer[pp][1], i1 = match_ouTer[pp][2], j1 = match_ouTer[pp][3], i2 = match_ouTer[pp][4], j2 = match_ouTer[pp][5], reversed = match_ouTer[pp][6], lonlat = match_ouTer[pp][7];

            for (int idx = 0; idx < NX; idx++) {
                int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                int I2_1 = i2 == -1 ? reversed ? checkIP_ouTer[NX-1-idx][0] : checkIP_ouTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? checkIP_ouTer[NY-1-idx][0] : checkIP_ouTer[idx][0] : j2;
                int I2_2 = i2 == -1 ? reversed ? checkIP_ouTer[NX-1-idx][1] : checkIP_ouTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? checkIP_ouTer[NY-1-idx][1] : checkIP_ouTer[idx][1] : j2;
                if (lonlat == 0) {
                    B = csswm[p1].lat[I1][J1];
                    A1 = csswm[p2].lat[I2_1][J2_1], A2 = csswm[p2].lat[I2_2][J2_2];
                }
                else {
                    B = csswm[p1].lon[I1][J1];
                    A1 = csswm[p2].lon[I2_1][J2_1], A2 = csswm[p2].lon[I2_2][J2_2];

                    if (A1 > A2 && (p1 == 0 || p2 == 0))  A2 += 2 * M_PI;
                    if (A1 > B && B < A2) B += 2 * M_PI;
                }

                alpha_tmp = alpha2D[I1][J1];
                beta_tmp = beta2D[I1][J1];

                alpha_B_tmp = interpolate(A1, A2, alpha2D[I2_1][J2_1], alpha2D[I2_2][J2_2], B);
                beta_B_tmp = interpolate(A1, A2, beta2D[I2_1][J2_1], beta2D[I2_2][J2_2], B);

                get_gLower(gLower_tmp, alpha_tmp, beta_tmp);
                get_IA(IA_tmp, p1, alpha_tmp, beta_tmp);
                get_A(A_tmp, p2, alpha_B_tmp, beta_B_tmp);
                get_gUpper(gUpper_tmp, alpha_B_tmp, beta_B_tmp);

                matrixMul(gLower_tmp, IA_tmp, tmp1);
                matrixMul(A_tmp, gUpper_tmp, tmp2);
                double tmp1_4[4], tmp2_4[4];
                count = 0;
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        tmp1_4[count] = tmp1[i][j];
                        tmp2_4[count] = tmp2[i][j];
                        count++;
                    }
                }
                matrixMul(tmp1_4, tmp2_4, mult);
                
                if (i1 == -1) {
                    if (j1 == 1) {
                        count = 0;
                        for (int i = 0; i < 2; i++) {
                            for (int j = 0; j < 2; j++) {
                                csswm[p1].IP_ouTer_D[I1][count] = mult[i][j];
                                count++;
                            }
                        }
                    }
                    else if (j1 == NY-2) {
                        count = 0;
                        for (int i = 0; i < 2; i++) {
                            for (int j = 0; j < 2; j++) {
                                csswm[p1].IP_ouTer_U[I1][count] = mult[i][j];
                                count++;
                            }
                        }
                    }
                }
                else if (j1 == -1) {
                    if (i1 == 1) {
                        count = 0;
                        for (int i = 0; i < 2; i++) {
                            for (int j = 0; j < 2; j++) {
                                csswm[p1].IP_ouTer_L[J1][count] = mult[i][j];
                                count++;
                            }
                        }
                    }
                    else if (i1 == NX-2) {
                        count = 0;
                        for (int i = 0; i < 2; i++) {
                            for (int j = 0; j < 2; j++) {
                                csswm[p1].IP_ouTer_R[J1][count] = mult[i][j];
                                count++;
                            }
                        }
                    }
                }
            }
        }
    #endif
}

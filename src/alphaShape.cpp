/******************************************************************************
 * Principal function:                                                         *
 *    Molecule shape overlay: Hard Sphere Model, Gaussian Sphere Model, flexible
 * 			     overlay
 *                                                                             *
 * Usage:
 *
 * Authors:                                                                    *
 *    James Li, SciNet Technologies, May 2012
 *    Xin Yan, May 2012
 *
 *                                                                             *
 ******************************************************************************/
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

#include <math.h>
#include <time.h>
#include "mathUtil.hpp"
#include "alphaShape.hpp"

static double gParam = 2.8284;

double overlay_hack(int NA1, double* xyz1, double* alpha1, double* w1, int NA2,
                                  double* xyz2, double* alpha2, double* w2, int nFeat1,
                                  int* FeatTypes1, double* featXYZ1, int nDirc1, int* dircTypes1,
                                  double* dircXYZ1, int nFeat2, int* FeatTypes2, double* featXYZ2,
                                  int nDirc2, int* dircTypes2, double* dircXYZ2, 
                                  bool equal_alpha, double* W, double* cVol, double* diffVolume); 
//
// Solve linear equation of dimension size 6 using Gaussian Elimination method
//
// hess*dx = der6
//
double LSolver6(double hess[6][6], double der6[6], double dx[6]) {
    double a[21], derv[6];
    double dmax;
    double base = 5.0;
    int nfc;
    int loop = 0;
    int ii[6];
    ii[0] = 0;
    for (int i = 1; i < 6; i++) ii[i] = ii[i - 1] + 7 - i;

    do {
        for (int i = 0; i < 6; i++) {
            derv[i] = der6[i];
            for (int j = i; j < 6; j++) a[ii[i] + j - i] = hess[i][j];
            a[ii[i]] += loop * base;
        }
        nfc = 0;
        for (int i = 0; i < 5; i++) {
            if (a[ii[i]] < 0) nfc++;
            //
            //   Do elimination
            //
            for (int j = i + 1; j < 6; j++) {
                double ex = a[ii[i] + j - i] / a[ii[i]];
                derv[j] -= derv[i] * ex;
                for (int k = j; k < 6; k++) a[ii[j] + k - j] -= a[ii[i] + k - i] * ex;
            }
        }
        if (a[20] < 0) nfc++;
        if (nfc > 0) {
            //   printf("Nagative eigenvalues !\n");
            //   printf("base = %lf !\n", base);
        }
        loop = 1;
        base *= 2.0;
    } while (nfc > 0);
    dmax = 0;
    for (int j = 5; j > -1; j--) {
        double tmp = derv[j];
        for (int k = j + 1; k < 6; k++) {
            tmp -= a[ii[j] + k - j] * dx[k];
        }
        dx[j] = tmp / a[ii[j]];
        if (fabs(dx[j]) > dmax) dmax = fabs(dx[j]);
    }
    return dmax;
}

void test_overlay()
{
     int NA1;
     double* xyz1;
     double* alpha1;
     double* w1;
     int NA2;
     double* xyz2;
     double* alpha2;
     double* w2;
     int nFeat1;
     int* FeatTypes1;
     double* featXYZ1;
     int nDirc1;
     int* dircTypes1;
     double* dircXYZ1;
     int nFeat2;
     int* FeatTypes2;
     double* featXYZ2;
     int nDirc2;
     int* dircTypes2;
     double* dircXYZ2;
     //bool trans;
     //bool equalR;
     double cVol[10];
     double diffVolume[10];
     double W[10]{1.0,
                  1.0,   // HBA in ownPharma
                  1.0,   // HBD
                  0.5,   // ARO
                  2.0,   // POS
                  2.0,   // NEG
                  2.0,   // HYD
                  1.0,   // for direction of HBA in overlay
                  1.0,   // for direction of HBD
                  0.5};  // for direction of ARO

    std::ifstream fin;
    fin.open("test_data.txt");
    std::string data_str;

    getline(fin, data_str);
    fin.close();
    std::istringstream fin0(data_str);

    fin0 >> NA1;
    alpha1 = new double[NA1];
    w1 = new double[NA1];
    xyz1 = new double[NA1*3];

    for (int i=0; i < NA1; i++) {
	int i3 = i*3;
	fin0 >> alpha1[i];
       	fin0 >> w1[i];
	fin0 >> xyz1[i3];
	fin0 >> xyz1[i3+1];
	fin0 >> xyz1[i3+2];
    }
    fin0 >> NA2;
    alpha2 = new double[NA2];
    w2 = new double[NA2];
    xyz2 = new double[NA2*3];
    for (int i=0; i < NA2; i++) {
	int i3 = i*3;
	fin0 >> alpha2[i];
       	fin0 >> w2[i];
	fin0 >> xyz2[i3];
	fin0 >> xyz2[i3+1];
	fin0 >> xyz2[i3+2];
    }
    fin0 >> nFeat1;
    FeatTypes1 = new int[nFeat1];
    featXYZ1 = new double[nFeat1*3];
    for (int i=0; i < nFeat1; i++) {
	int i3 = i*3;
	fin0 >> FeatTypes1[i];
	fin0 >> featXYZ1[i3];
	fin0 >> featXYZ1[i3+1];
	fin0 >> featXYZ1[i3+2];
    }
    fin0 >> nDirc1;
    dircTypes1 = new int[nDirc1];
    dircXYZ1 = new double[nDirc1*3];
    for (int i=0; i < nDirc1; i++) {
	int i3 = i*3;
	fin0 >> dircTypes1[i];
	fin0 >> dircXYZ1[i3];
	fin0 >> dircXYZ1[i3+1];
	fin0 >> dircXYZ1[i3+2];
    }
    fin0 >> nFeat2;
    FeatTypes2 = new int[nFeat2];
    featXYZ2 = new double[nFeat2*3];
    for (int i=0; i < nFeat2; i++) {
	int i3 = i*3;
	fin0 >> FeatTypes2[i];
	fin0 >> featXYZ2[i3];
	fin0 >> featXYZ2[i3+1];
	fin0 >> featXYZ2[i3+2];
    }
    fin0 >> nDirc2;
    dircTypes2 = new int[nDirc2];
    dircXYZ2 = new double[nDirc2*3];
    for (int i=0; i < nDirc2; i++) {
	int i3 = i*3;
	fin0 >> dircTypes2[i];
	fin0 >> dircXYZ2[i3];
	fin0 >> dircXYZ2[i3+1];
	fin0 >> dircXYZ2[i3+2];
    }

    bool equal_alpha = false;

    double* xyzA = new double[NA1*3];
    double* xyzB = new double[NA2*3];
    double* featXYZA = new double[nFeat1*3];
    double* featXYZB = new double[nFeat2*3];
    double* dircXYZA = new double[nDirc1*3];
    double* dircXYZB = new double[nDirc2*3];
    printf("NA1 = %d, NA2 = %d\n",NA1, NA2);
    for (int i=0; i < 10000; i++) {
        for (int k=0; k < NA1*3; k++) {
            xyzA[k] = xyz1[k];
        }
        for (int k=0; k < NA2*3; k++) {
            xyzB[k] = xyz2[k];
        }
        for (int k=0; k < nFeat1*3; k++) {
            featXYZA[k] = featXYZ1[k];
        }
        for (int k=0; k < nFeat2*3; k++) {
            featXYZB[k] = featXYZ2[k];
        }
        for (int k=0; k < nDirc1*3; k++) {
            dircXYZA[k] = dircXYZ1[k];
        }
        for (int k=0; k < nDirc2*3; k++) {
            dircXYZB[k] = dircXYZ2[k];
        }
        overlay_hack(NA1, xyzA, alpha1, w1, 
                     NA2, xyzB, alpha2, w2, 
                     nFeat1, FeatTypes1, featXYZA, 
                     nDirc1, dircTypes1, dircXYZA, 
                     nFeat2, FeatTypes2, featXYZB, 
                     nDirc2, dircTypes2, dircXYZB,
		             equal_alpha, W, cVol, diffVolume);
   }
    
};

double overlay_hack(int NA1, double* xyz1, double* alpha1, double* w1, int NA2,
                                  double* xyz2, double* alpha2, double* w2, int nFeat1,
                                  int* FeatTypes1, double* featXYZ1, int nDirc1, int* dircTypes1,
                                  double* dircXYZ1, int nFeat2, int* FeatTypes2, double* featXYZ2,
                                  int nDirc2, int* dircTypes2, double* dircXYZ2,
                                  bool equal_alpha, double* W, double* cVol, double* diffVolume) 
{
    //
    // r' = R*r + T
    //
    // r': new position
    // r:  Original position of mol1
    // R:  Rotation matrix
    // T:  Translation
    //
    int ij;
    int i;
    int i3,i6;
    double pi = 3.14159265358979;
    double p2 = gParam * gParam;
    double dltn;
    double pidltn;

    double molvol;
    double molvol2;
    double diffVol;
    double xyz0[3], xyzr[3];
    double r[3][3];
    double xyz[3000];
    double T[3];
    double a, b, g;
    //
    // a, b, g, x, y, z
    //
    double hess[6][6], dtr[6], dx[6];
    double Thresh = 0.01;
    int maxIter = 25;
    double dmax = 0;
    double hessOld[6][6];

    double* cderv = new double[NA1 * 3];
    double* chess = new double[NA1 * 6];
    // Save old coordinates
    // Set zeros for cderv and chess
    for (i = 0; i < NA1; i++) {
        i3 = i * 3;
        i6 = i3 + i3;
        xyz[i3] = xyz1[i3];
        xyz[i3 + 1] = xyz1[i3 + 1];
        xyz[i3 + 2] = xyz1[i3 + 2];
    }

    // Apply transformation
    // Initial guess

    T[0] = 0;
    T[1] = 0;
    T[2] = 0;
    a = 0;
    b = 0;
    g = 0;

    int iter = 0;
    molvol2 = 0;
    double damp = 0.5;

    double* wxpp = new double[NA1 * NA2];
    double* afac = new double[NA1 * NA2];
    double* ahaf = new double[NA1 * NA2];

    if (!equal_alpha)  // all vdwR are equal, so all alphas are equal
    {
        for (i = 0; i < NA1; i++) {
            double w1p2 = w1[i] * p2;
            int ina2 = i * NA2;
            for (int j = 0; j < NA2; j++) {
                dltn = alpha1[i] + alpha2[j];
                pidltn = pi / dltn;
                wxpp[ina2 + j] = w1p2 * w2[j] * pidltn * sqrt(pidltn);
                afac[ina2 + j] = alpha1[i]/(alpha1[i]+alpha2[j]);
                ahaf[ina2 + j] = afac[ina2 + j]*alpha2[j];
            }
        }
    }

    do {
        iter++;
        for (i = 0; i < 6; i++) {
            dtr[i] = 0;
            for (int j = 0; j < 6; j++) {
                hess[i][j] = 0;
            }
        }
        for (i = 0; i < NA1; i++) {
            i3 = i * 3;
            i6 = i * 6;
            cderv[i3] = 0;
            cderv[i3+1] = 0;
            cderv[i3+2] = 0;
            chess[i6] = 0;
            chess[i6+1] = 0;
            chess[i6+2] = 0;
            chess[i6+3] = 0;
            chess[i6+4] = 0;
            chess[i6+5] = 0;
        }
        std::cout << "Iter = " << iter << "\n";
        if (equal_alpha)
            molvol = cal_cart_dervs(NA1, xyz1, w1, NA2, xyz2, w2, cderv, chess);
        else
            molvol = cal_cart_dervs2(NA1, xyz1, alpha1, w1, NA2, xyz2, alpha2, w2, 
                                        cderv, chess, wxpp, afac, ahaf);         
        cart_dervs_2rt(NA1, xyz1, cderv, chess, dtr, hess);

        cVol[0] = molvol;
        molvol = molvol * W[0];
        for (i = 0; i < 6; i++) {
            dtr[i] = dtr[i] * W[0];
            for (int j = 0; j < 6; j++) {
                hess[i][j] = hess[i][j] * W[0];
            }
        }

        for (i = 1; i < 10; i++) cVol[i] = 0;

        int feat_num1[10];
        int feat_num2[10];
        for (int i=0; i < 10; i++) {
            feat_num1[i] = 0;
            feat_num2[i] = 0;
        }
        for (int i=0; i < nFeat1; i++) {
            feat_num1[FeatTypes1[i]]++;
        }
        for (int i=0; i < nFeat2; i++) {
            feat_num2[FeatTypes2[i]]++;
        }   
        for (int i=0; i < nDirc1; i++) {
            feat_num1[6+dircTypes1[i]]++;
        }
        for (int i=0; i < nDirc2; i++) {
            feat_num2[6+dircTypes2[i]]++;
        }
        for (int ftype=1; ftype < 10; ftype++) {
            int na1 = feat_num1[ftype];
            int na2 = feat_num2[ftype];
            int nSize = na1*na2;
            if (nSize == 0) continue;
            double f_xyz1[na1*3];
            double f_xyz2[na2*3];
            int nxyz3 = 0;
            for (int i=0; i < nFeat1; i++) {
                if (FeatTypes1[i] == ftype) {
                    i3 = i*3;
                    f_xyz1[nxyz3] = featXYZ1[i3];
                    f_xyz1[nxyz3+1] = featXYZ1[i3+1];
                    f_xyz1[nxyz3+2] = featXYZ1[i3+2];
                    nxyz3 += 3;
                } 
            } 
            nxyz3 = 0;
            for (int i=0; i < nFeat2; i++) {
                if (FeatTypes2[i] == ftype) {
                    i3 = i*3;
                    f_xyz2[nxyz3] = featXYZ2[i3];
                    f_xyz2[nxyz3+1] = featXYZ2[i3+1];
                    f_xyz2[nxyz3+2] = featXYZ2[i3+2];
                    nxyz3 += 3;
                }
            }
            nxyz3 = 0;
            for (int i=0; i < nDirc1; i++) {
                if (dircTypes1[i] == ftype-6) {
                    i3 = i*3;
                    f_xyz1[nxyz3] = dircXYZ1[i3];
                    f_xyz1[nxyz3+1] = dircXYZ1[i3+1];
                    f_xyz1[nxyz3+2] = dircXYZ1[i3+2];
                    nxyz3 += 3;
                } 
            } 
            nxyz3 = 0;
            for (int i=0; i < nDirc2; i++) {
                if (dircTypes2[i] == ftype-6) {
                    i3 = i*3;
                    f_xyz2[nxyz3] = dircXYZ2[i3];
                    f_xyz2[nxyz3+1] = dircXYZ2[i3+1];
                    f_xyz2[nxyz3+2] = dircXYZ2[i3+2];
                    nxyz3 += 3;
                }
            }
            double wf1[na1];
            double wf2[na2];
            double cderv1[na1*3];
            double chess1[na1*6];
            for (int i=0; i < na1; i++) {
                wf1[i] = W[ftype];
            }
            for (int i=0; i < na2; i++) {
                wf2[i] = 1.0;
            }
            for (i = 0; i < na1; i++) {
                i3 = i * 3;
                i6 = i * 6;
                cderv1[i3] = 0;
                cderv1[i3+1] = 0;
                cderv1[i3+2] = 0;
                chess1[i6] = 0;
                chess1[i6+1] = 0;
                chess1[i6+2] = 0;
                chess1[i6+3] = 0;
                chess1[i6+4] = 0;
                chess1[i6+5] = 0;
            }
            cVol[ftype] = cal_cart_dervs(na1, f_xyz1, wf1, na2, f_xyz2, wf2, cderv1, chess1);
            molvol += cVol[ftype];
            cart_dervs_2rt(na1, f_xyz1, cderv1, chess1, dtr, hess);
        }
        double dmaxFactor = 1.0;
        diffVol = fabs(molvol - molvol2);
        dmax = LSolver6(hess, dtr, dx);

        damp = 1.0;
        if (molvol < molvol2) {
            dmaxFactor = 0.5;
        }

        if (dmax > 0.5 * dmaxFactor) damp = 0.5 * dmaxFactor / dmax;
        molvol2 = molvol;
        a = -damp * dx[0];
        b = -damp * dx[1];
        g = -damp * dx[2];
        T[0] = -damp * dx[3];
        T[1] = -damp * dx[4];
        T[2] = -damp * dx[5];
        // Translation
        RotMat(a, b, g, r);
        for (i = 0; i < NA1; i++) {
            i3 = i * 3;
            Rotation(xyz1 + i3, r, xyzr);
            xyz1[i3] = xyzr[0] + T[0];
            xyz1[i3 + 1] = xyzr[1] + T[1];
            xyz1[i3 + 2] = xyzr[2] + T[2];
        }

        for (i = 0; i < nFeat1; i++) {
            i3 = i * 3;
            Rotation(featXYZ1 + i3, r, xyzr);
            featXYZ1[i3] = xyzr[0] + T[0];
            featXYZ1[i3 + 1] = xyzr[1] + T[1];
            featXYZ1[i3 + 2] = xyzr[2] + T[2];
        }

        for (i = 0; i < nDirc1; i++) {
            i3 = i * 3;
            Rotation(dircXYZ1 + i3, r, xyzr);
            dircXYZ1[i3] = xyzr[0] + T[0];
            dircXYZ1[i3 + 1] = xyzr[1] + T[1];
            dircXYZ1[i3 + 2] = xyzr[2] + T[2];
        }

    } while (diffVol > Thresh && iter < maxIter);
    // printf("diffVol: %f; Iter: %d\n",diffVol,iter);
    if (0 != diffVolume) diffVolume[0] = diffVol;

    double molvx = cVol[0];
    for (i = 1; i < 10; i++) {
        if (0.001 > W[i]) continue;
        cVol[i] = cVol[i] / W[i];
	molvx += cVol[i];
    }

    // Copy xyz back to mol
    for (i = 0; i < NA1; i++) {
        i3 = i * 3;
        xyz1[i3] = xyz[i3];
        xyz1[i3 + 1] = xyz[i3 + 1];
        xyz1[i3 + 2] = xyz[i3 + 2];
    }

    delete[] cderv;
    delete[] chess;
    delete[] wxpp;
    std::cout << "molvol, molvx = " << molvol << " " << molvx << "\n";
    return molvol;
}

double cal_cart_dervs(int na1, double* xyz1, double* w1, int na2, double* xyz2,
                         double* w2, double* cderv, double* chess) {
    int i, j, ina2;
    int i3, i31, i32, i6, j3, j31, j32;
    double kn, dist;
    double tmpd, tav, tavad, tavap;
    double crossVol = 0;
    double rn[3];
    double derv[3], hess[6];    
    double alpha = 0.74628778;  // gFactor/(VDWR*VDWR);
    double alphaH = alpha / 2;
    double alphaD = alpha * 2;

    for (i = 0; i < na1; i++) {
        i3 = i * 3;
        i31 = i3+1;
        i32 = i3+2;
        i6 = i3 + i3;
        ina2 = i * na2;
        double wi = 24.4287905*w1[i];
        derv[0] = 0;
        derv[1] = 0;
        derv[2] = 0;
        hess[0] = 0;
        hess[1] = 0;
        hess[2] = 0;
        hess[3] = 0;
        hess[4] = 0;
        hess[5] = 0;

        for (j = 0; j < na2; j++) {
            j3 = j * 3;
            j31 = j3+1;
            j32 = j3+2;
            dist = (xyz1[i3] - xyz2[j3]) * (xyz1[i3] - xyz2[j3]) +
                    (xyz1[i31] - xyz2[j31]) * (xyz1[i31] - xyz2[j31]) +
                    (xyz1[i32] - xyz2[j32]) * (xyz1[i32] - xyz2[j32]);
            rn[0] = (xyz1[i3] + xyz2[j3]) * 0.5;
            rn[1] = (xyz1[i31] + xyz2[j31]) * 0.5;
            rn[2] = (xyz1[i32] + xyz2[j32]) * 0.5;
            kn = exp(-alphaH * dist);
            tmpd =  wi*w2[j]*kn;
            crossVol += tmpd;
            // First derivatives
            tav = alphaD * tmpd;

            double dv0 = xyz1[i3] - rn[0];
            double dv1 = xyz1[i3+1] - rn[1];
            double dv2 = xyz1[i3+2] - rn[2];

            derv[0] += tav * dv0;
            derv[1] += tav * dv1;
            derv[2] += tav * dv2;
            tavad = tav * (-0.5);
            tavap = (tav + tav) * alpha;
            // Hessian
            hess[0] -= tavad + tavap * dv0 * dv0;
            hess[1] -= tavap * dv0 * dv1;
            hess[2] -= tavad + tavap * dv1 * dv1;
            hess[3] -= tavap * dv0 * dv2;
            hess[4] -= tavap * dv1 * dv2;
            hess[5] -= tavad + tavap * dv2 * dv2;
        }
        // save dervs of each atom in na1
        cderv[i3]     += derv[0];
        cderv[i3 + 1] += derv[1];
        cderv[i3 + 2] += derv[2];
        chess[i6]     += hess[0];
        chess[i6 + 1] += hess[1];
        chess[i6 + 2] += hess[2];
        chess[i6 + 3] += hess[3];
        chess[i6 + 4] += hess[4];
        chess[i6 + 5] += hess[5];
    }
    return crossVol;
}

double cal_cart_dervs2(int na1, double* xyz1, double* alpha1, double* w1, int na2, double* xyz2,
                         double* alpha2, double* w2, double* cderv, double* chess, 
                         double* wxpp, double* afac, double* ahaf) {
    int i, j, ina2;
    int i3, i31, i32, i6, j3, j31, j32; 
    double kn, dist;
    double tmpd, tav, tavad, tavap;
    double crossVol = 0;
    double rn[3];
    double derv[3], hess[6];
    for (i = 0; i < na1; i++) {
        i3 = i * 3;
        i31 = i3+1;
        i32 = i3+2;
        i6 = i3 + i3;
        ina2 = i * na2;
        derv[0] = 0;
        derv[1] = 0;
        derv[2] = 0;
        hess[0] = 0;
        hess[1] = 0;
        hess[2] = 0;
        hess[3] = 0;
        hess[4] = 0;
        hess[5] = 0;
        double a12 = alpha1[i]*2;
        for (j = 0; j < na2; j++) {
            j3 = j * 3;
            j31 = j3+1;
            j32 = j3+2;
            dist = (xyz1[i3] - xyz2[j3]) * (xyz1[i3] - xyz2[j3]) +
                    (xyz1[i31] - xyz2[j31]) * (xyz1[i31] - xyz2[j31]) +
                    (xyz1[i32] - xyz2[j32]) * (xyz1[i32] - xyz2[j32]);
            rn[0] = xyz1[i3]*afac[ina2+j] + xyz2[j3]*(1.0-afac[ina2+j]);
            rn[1] = xyz1[i31]*afac[ina2+j] + xyz2[j3+1]*(1.0-afac[ina2+j]);
            rn[2] = xyz1[i32]*afac[ina2+j] + xyz2[j3+2]*(1.0-afac[ina2+j]);
            kn = exp(-ahaf[ina2+j] * dist);
            tmpd = wxpp[ina2 + j] * kn;
            crossVol += tmpd;
            tav = a12 * tmpd;
            double dv0 = xyz1[i3] - rn[0];
            double dv1 = xyz1[i31] - rn[1];
            double dv2 = xyz1[i32] - rn[2];

            // First derivatives
            derv[0] += tav * dv0;
            derv[1] += tav * dv1;
            derv[2] += tav * dv2;
            tavad = tav * (afac[ina2+j] - 1);
            tavap = tav * a12;

            // Hessian
            double tavap_dv0 = tavap*dv0;
            double tavap_dv1 = tavap*dv1;
            double tavap_dv2 = tavap*dv2;
            hess[0] -= tavad + tavap_dv0 * dv0;
            hess[1] -= tavap_dv0 * dv1;
            hess[2] -= tavad + tavap_dv1 * dv1;
            hess[3] -= tavap_dv0 * dv2;
            hess[4] -= tavap_dv1 * dv2;
            hess[5] -= tavad + tavap_dv2 * dv2;
        }
        // save dervs of each atom in na1
        cderv[i3]     += derv[0];
        cderv[i31]    += derv[1];
        cderv[i32]    += derv[2];
        chess[i6]     += hess[0];
        chess[i6 + 1] += hess[1];
        chess[i6 + 2] += hess[2];
        chess[i6 + 3] += hess[3];
        chess[i6 + 4] += hess[4];
        chess[i6 + 5] += hess[5];
    }
    return crossVol;
}

void cart_dervs_2rt(int na1, double* xyz1, double* cderv, double* chess, double rt_derv[6],
                       double rt_hess[6][6]) {
    int i, i3, i6;
    double x, y, z;
    double derv[3];
    double hess[6];
    double xx, xy, xz, yz, yy, zz;
    double derv0x, derv0y, derv0z;
    double derv1x, derv1y, derv1z;
    double derv2x, derv2y, derv2z;

    for (i = 0; i < na1; i++) {
        i3 = i * 3;
        i6 = i3 + i3;
        derv[0] = cderv[i3];
        derv[1] = cderv[i3 + 1];
        derv[2] = cderv[i3 + 2];
        hess[0] = chess[i6];
        hess[1] = chess[i6 + 1];
        hess[2] = chess[i6 + 2];
        hess[3] = chess[i6 + 3];
        hess[4] = chess[i6 + 4];
        hess[5] = chess[i6 + 5];

        x = xyz1[i3];
        y = xyz1[i3 + 1];
        z = xyz1[i3 + 2];
        xx = x*x;
        xy = x*y;
        xz = x*z;
        yz = y*z;
        yy = y*y;
        zz = z*z;

        derv0x = derv[0]*x;
        derv0y = derv[0]*y;
        derv0z = derv[0]*z;
        derv1x = derv[1]*x;
        derv1y = derv[1]*y;
        derv1z = derv[1]*z;
        derv2x = derv[2]*x;
        derv2y = derv[2]*y;
        derv2z = derv[2]*z;
        //
        // First derivatives
        //
        rt_derv[3] += derv[0];
        rt_derv[4] += derv[1];
        rt_derv[5] += derv[2];
        rt_derv[0] += -derv0y + derv1x;
        rt_derv[1] += -derv1z + derv2y;
        rt_derv[2] +=  derv0z - derv2x;
        //
        // Hessian. variables 0=alpha, 1=beta, 2=gamma, 3=tx, 4=ty, 5=tx
        //
        rt_hess[3][3] += hess[0];
        rt_hess[3][4] += hess[1];
        rt_hess[4][4] += hess[2];
        rt_hess[3][5] += hess[3];
        rt_hess[4][5] += hess[4];
        rt_hess[5][5] += hess[5];
        rt_hess[0][0] += -derv0x - derv1y +
                         hess[0] * yy -
                         hess[1] * xy * 2.0 + hess[2] * xx;
        rt_hess[0][1] +=  derv0z + hess[1] * yz -
                         hess[2] * xz - hess[3] * yy +
                         hess[4] * xy;
        rt_hess[1][1] += -derv1y - derv2z +
                         hess[2] * zz -
                         hess[4] * yz * 2.0 +
                         hess[5] * yy;
        rt_hess[0][2] +=  derv2y - hess[0] * yz +
                         hess[1] * xz + hess[3] * xy -
                         hess[4] * xx;
        rt_hess[1][2] +=  derv0y - hess[1] * zz +
                         hess[3] * yz + hess[4] * xz -
                         hess[5] * xy;
        rt_hess[2][2] += -derv0x - derv2z +
                         hess[0] * zz -
                         hess[3] * xz * 2.0 + hess[5] * xx;
        rt_hess[0][3] += -hess[0] * y + hess[1] * x;
        rt_hess[1][3] += -hess[1] * z + hess[3] * y;
        rt_hess[2][3] +=  hess[0] * z - hess[3] * x;

        rt_hess[0][4] += -hess[1] * y + hess[2] * x;
        rt_hess[1][4] += -hess[2] * z + hess[4] * y;
        rt_hess[2][4] +=  hess[1] * z - hess[4] * x;

        rt_hess[0][5] += -hess[3] * y + hess[4] * x;
        rt_hess[1][5] += -hess[4] * z + hess[5] * y;
        rt_hess[2][5] +=  hess[3] * z - hess[5] * x;
    }
    for (i = 0; i < 6; i++) {
        for (i3 = 0; i3 < i + 1; i3++) {
            rt_hess[i][i3] = rt_hess[i3][i];
        }
    }
}
int main()
{
    test_overlay();
    return 0;
}

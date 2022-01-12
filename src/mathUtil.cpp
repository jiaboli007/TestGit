/******************************************************************************
 *    Description: Implementation of Molecular Foundation Classes              *
 *                                                                             *
 *    Author:      James Li                                                    *
 *                                                                             *
 *    Date:        March 2011                                                  *
 *                                                                             *
 ******************************************************************************/
#include <cmath>
#include <iostream>
#include "mathUtil.hpp"

void Rotation(double xyz0[3], double r[3][3], double xyz1[3]) {
    xyz1[0] = r[0][0] * xyz0[0] + r[0][1] * xyz0[1] + r[0][2] * xyz0[2];
    xyz1[1] = r[1][0] * xyz0[0] + r[1][1] * xyz0[1] + r[1][2] * xyz0[2];
    xyz1[2] = r[2][0] * xyz0[0] + r[2][1] * xyz0[1] + r[2][2] * xyz0[2];
}

void Rotation_float(float xyz0[3], float r[3][3], float xyz1[3]) {
    xyz1[0] = r[0][0] * xyz0[0] + r[0][1] * xyz0[1] + r[0][2] * xyz0[2];
    xyz1[1] = r[1][0] * xyz0[0] + r[1][1] * xyz0[1] + r[1][2] * xyz0[2];
    xyz1[2] = r[2][0] * xyz0[0] + r[2][1] * xyz0[1] + r[2][2] * xyz0[2];
}

//
//
//   Rotation matrix
//   y(gama)*z(alpha)*x(beta)
//
//   degree should be radian
//
void RotMat(double alpha, double beta, double gama, double r[3][3]) {
    double sa = sin(alpha);
    double ca = cos(alpha);
    double sb = sin(beta);
    double cb = cos(beta);
    double sh = sin(gama);
    double ch = cos(gama);
    // Rotation matrix
    r[0][0] = ch * ca;
    r[1][0] = sa;
    r[2][0] = -sh * ca;
    r[0][1] = -ch * sa * cb + sh * sb;
    r[1][1] = ca * cb;
    r[2][1] = sh * sa * cb + ch * sb;
    r[0][2] = ch * sa * sb + sh * cb;
    r[1][2] = -ca * sb;
    r[2][2] = -sh * sa * sb + ch * cb;
}

void RotMat_float(float alpha, float beta, float gama, float r[3][3]) {
    float sa = sin(alpha);
    float ca = cos(alpha);
    float sb = sin(beta);
    float cb = cos(beta);
    float sh = sin(gama);
    float ch = cos(gama);
    // Rotation matrix
    r[0][0] = ch * ca;
    r[1][0] = sa;
    r[2][0] = -sh * ca;
    r[0][1] = -ch * sa * cb + sh * sb;
    r[1][1] = ca * cb;
    r[2][1] = sh * sa * cb + ch * sb;
    r[0][2] = ch * sa * sb + sh * cb;
    r[1][2] = -ca * sb;
    r[2][2] = -sh * sa * sb + ch * cb;
}

//
//
//   Rotation matrix
//   y(gama)*z(alpha)*x(beta)
//
//   Return: the rational matrix, the first derivative matrices,
//           and the second derivative matrices
void Euler2Mat(double alpha, double beta, double gama, double r[3][3], double DRA[3][3],
               double DRB[3][3], double DRG[3][3], double DDRAB[3][3], double DDRAG[3][3],
               double DDRBG[3][3], double DDRAA[3][3], double DDRBB[3][3], double DDRGG[3][3]) {
    double sa = sin(alpha);
    double ca = cos(alpha);
    double sb = sin(beta);
    double cb = cos(beta);
    double sh = sin(gama);
    double ch = cos(gama);
    // Rotation matrix
    r[0][0] = ch * ca;
    r[1][0] = sa;
    r[2][0] = -sh * ca;
    r[0][1] = -ch * sa * cb + sh * sb;
    r[1][1] = ca * cb;
    r[2][1] = sh * sa * cb + ch * sb;
    r[0][2] = ch * sa * sb + sh * cb;
    r[1][2] = -ca * sb;
    r[2][2] = -sh * sa * sb + ch * cb;
    // alpha
    DRA[0][0] = -ch * sa;
    DRA[1][0] = ca;
    DRA[2][0] = sh * sa;
    DRA[0][1] = -ch * ca * cb;
    DRA[1][1] = -sa * cb;
    DRA[2][1] = sh * ca * cb;
    DRA[0][2] = ch * ca * sb;
    DRA[1][2] = sa * sb;
    DRA[2][2] = -sh * ca * sb;
    // beta
    DRB[0][0] = 0;
    DRB[1][0] = 0;
    DRB[2][0] = 0;
    DRB[0][1] = ch * sa * sb + sh * cb;
    DRB[1][1] = -ca * sb;
    DRB[2][1] = -sh * sa * sb + ch * cb;
    DRB[0][2] = ch * sa * cb - sh * sb;
    DRB[1][2] = -ca * cb;
    DRB[2][2] = -sh * sa * cb - ch * sb;
    // gama
    DRG[0][0] = -sh * ca;
    DRG[1][0] = 0;
    DRG[2][0] = -ch * ca;
    DRG[0][1] = sh * sa * cb + ch * sb;
    DRG[1][1] = 0;
    DRG[2][1] = ch * sa * cb - sh * sb;
    DRG[0][2] = -sh * sa * sb + ch * cb;
    DRG[1][2] = 0;
    DRG[2][2] = -ch * sa * sb - sh * cb;
    // alpah-beta
    DDRAB[0][0] = 0;
    DDRAB[1][0] = 0;
    DDRAB[2][0] = 0;
    DDRAB[0][1] = ch * ca * sb;
    DDRAB[1][1] = sa * sb;
    DDRAB[2][1] = -sh * ca * sb;
    DDRAB[0][2] = ch * ca * cb;
    DDRAB[1][2] = sa * cb;
    DDRAB[2][2] = -sh * ca * cb;
    // alpha-gama
    DDRAG[0][0] = sh * sa;
    DDRAG[1][0] = 0;
    DDRAG[2][0] = ch * sa;
    DDRAG[0][1] = sh * ca * cb;
    DDRAG[1][1] = 0;
    DDRAG[2][1] = ch * ca * cb;
    DDRAG[0][2] = -sh * ca * sb;
    DDRAG[1][2] = 0;
    DDRAG[2][2] = -ch * ca * sb;
    // beat-gama
    DDRBG[0][0] = 0;
    DDRBG[1][0] = 0;
    DDRBG[2][0] = 0;
    DDRBG[0][1] = -sh * sa * sb + ch * cb;
    DDRBG[1][1] = 0;
    DDRBG[2][1] = -ch * sa * sb - sh * cb;
    DDRBG[0][2] = -sh * sa * cb - ch * sb;
    DDRBG[1][2] = 0;
    DDRBG[2][2] = -ch * sa * cb + sh * sb;
    // alpha-alpha

    DDRAA[0][0] = -ch * ca;
    DDRAA[1][0] = -sa;
    DDRAA[2][0] = sh * ca;
    DDRAA[0][1] = ch * sa * cb;
    DDRAA[1][1] = -ca * cb;
    DDRAA[2][1] = -sh * sa * cb;
    DDRAA[0][2] = -ch * sa * sb;
    DDRAA[1][2] = ca * sb;
    DDRAA[2][2] = sh * sa * sb;
    // beta-beta
    DDRBB[0][0] = 0;
    DDRBB[1][0] = 0;
    DDRBB[2][0] = 0;
    DDRBB[0][1] = ch * sa * cb - sh * sb;
    DDRBB[1][1] = -ca * cb;
    DDRBB[2][1] = -sh * sa * cb - ch * sb;
    DDRBB[0][2] = -ch * sa * sb - sh * cb;
    DDRBB[1][2] = ca * sb;
    DDRBB[2][2] = sh * sa * sb - ch * cb;
    // gama-gama
    DRG[0][0] = -sh * ca;
    DRG[1][0] = 0;
    DRG[2][0] = -ch * ca;
    DRG[0][1] = sh * sa * cb + ch * sb;
    DRG[1][1] = 0;
    DRG[2][1] = ch * sa * cb - sh * sb;
    DRG[0][2] = -sh * sa * sb + ch * cb;
    DRG[1][2] = 0;
    DRG[2][2] = -ch * sa * sb - sh * cb;

    DDRGG[0][0] = -ch * ca;
    DDRGG[1][0] = 0;
    DDRGG[2][0] = sh * ca;
    DDRGG[0][1] = ch * sa * cb - sh * sb;
    DDRGG[1][1] = 0;
    DDRGG[2][1] = -sh * sa * cb - ch * sb;
    DDRGG[0][2] = -ch * sa * sb - sh * cb;
    DDRGG[1][2] = 0;
    DDRGG[2][2] = sh * sa * sb - ch * cb;
}

Quaternion::Quaternion() {
    x = 0.0;
    y = 0.0;
    z = 0.0;
    w = 1.0;
}

Quaternion::Quaternion(Quaternion* q) {
    x = q->x;
    y = q->y;
    z = q->z;
    w = q->w;
}

//
// angle should be radian
//
Quaternion::Quaternion(double angle, double v[3]) {
    double cosa = cos(angle / 2.0);
    double sina = sin(angle / 2.0);
    double factor = sina / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    x = factor * v[0];
    y = factor * v[1];
    z = factor * v[2];
    w = cosa;
}

void QuatToRMatrix(double m[3][3], Quaternion* Q) {
    double q[4];
    q[0] = Q->x;
    q[1] = Q->y;
    q[2] = Q->z;
    q[3] = Q->w;
    QuatToRMatrix(m, q);
}

void RMatrixToQuat(double m[3][3], Quaternion* Q) {
    double S = 0;
    double T = 1 + m[0][0] + m[1][1] + m[2][2];
    if (T > 0.00000001) {
        S = sqrt(T) * 2;
        Q->x = (m[2][1] - m[1][2]) / S;
        Q->y = (m[0][2] - m[2][0]) / S;
        Q->z = (m[1][0] - m[0][1]) / S;
        Q->w = S * 0.25;
    } else if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {
        S = sqrt(1.0 + m[0][0] - m[1][1] - m[2][2]) * 2;
        Q->x = S * 0.25;
        Q->y = (m[1][0] + m[0][1]) / S;
        Q->z = (m[0][2] + m[2][0]) / S;
        Q->w = (m[2][1] - m[1][2]) / S;
    } else if (m[1][1] > m[2][2]) {
        S = sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]) * 2;
        Q->x = (m[1][0] + m[0][1]) / S;
        Q->y = S * 0.25;
        Q->z = (m[2][1] + m[1][2]) / S;
        Q->w = (m[0][2] - m[2][0]) / S;
    } else {
        S = sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]) * 2;
        Q->x = (m[0][2] + m[2][0]) / S;
        Q->y = (m[2][1] + m[1][2]) / S;
        Q->z = S * 0.25;
        Q->w = (m[1][0] - m[0][1]) / S;
    }
}

void QuatToRMatrix(double m[3][3], double q[4]) {
    m[0][0] = 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);
    m[0][1] = 2.0 * (q[0] * q[1] - q[2] * q[3]);
    m[0][2] = 2.0 * (q[2] * q[0] + q[1] * q[3]);

    m[1][0] = 2.0 * (q[0] * q[1] + q[2] * q[3]);
    m[1][1] = 1.0 - 2.0 * (q[2] * q[2] + q[0] * q[0]);
    m[1][2] = 2.0 * (q[1] * q[2] - q[0] * q[3]);

    m[2][0] = 2.0 * (q[2] * q[0] - q[1] * q[3]);
    m[2][1] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
    m[2][2] = 1.0 - 2.0 * (q[1] * q[1] + q[0] * q[0]);
}

void QuatMultiple(Quaternion* Q1, Quaternion* Q2, Quaternion* Q3) {
    double q1[4], q2[4], q3[4];
    q1[0] = Q1->x;
    q1[1] = Q1->y;
    q1[2] = Q1->z;
    q1[3] = Q1->w;
    q2[0] = Q2->x;
    q2[1] = Q2->y;
    q2[2] = Q2->z;
    q2[3] = Q2->w;
    QuatMultiple(q1, q2, q3);
    Q3->x = q3[0];
    Q3->y = q3[1];
    Q3->z = q3[2];
    Q3->w = q3[3];
}

void QuatMultiple(double p1[4], double p2[4], double dest[4]) {
    double q10, q11, q12, q13, q20, q21, q22, q23;

    q11 = p1[0];
    q12 = p1[1];
    q13 = p1[2];
    q10 = p1[3];
    q21 = p2[0];
    q22 = p2[1];
    q23 = p2[2];
    q20 = p2[3];

    dest[3] = q10 * q20 - q11 * q21 - q12 * q22 - q13 * q23;
    dest[0] = q10 * q21 + q11 * q20 + q12 * q23 - q13 * q22;
    dest[1] = q10 * q22 + q12 * q20 + q13 * q21 - q11 * q23;
    dest[2] = q10 * q23 + q13 * q20 + q11 * q22 - q12 * q21;
}

bool QuatEqual(Quaternion* Q1, Quaternion* Q2) {
    double thres = 0.00000001;
    double dx = fabs(Q1->x - Q2->x);
    double dy = fabs(Q1->y - Q2->y);
    double dz = fabs(Q1->z - Q2->z);
    double dw = fabs(Q1->w - Q2->w);
    if (dx < thres && dy < thres && dz < thres && dw < thres) return true;

    return false;
}

int Jacobi(double* A, double* U, double* E, int NDIM, int N) {
    /*********************************************************************
     *   ARGUMENTS.                                                       *
     *        A       MATRIX TO BE DIAGONALIZED.                          *
     *        U       MATRIX OF EIGENVECTORS.                             *
     *        E       LIST OF EIGENVALUES                                 *
     *       NV       # OF EIGENVALUES NEEDED                             *
     *     NDIM       DIMENSION OF THE MATRICES AND THE VECTOR.           *
     *   Return       0: GOOD, 0: FAIL
     **********************************************************************/
    int MAXIT = 100, ITCOUNT = 0, JUP, I, J, K, II, JJ, KI, KJ, IJ, JI, IK, JK;
    double ZER = 0, ONE = 1, TWO = 2, FOR = 4, EPS = 1.0E-30;
    double AMAX, AII, AJJ, AOD, ASQ, DIFFR, SIGN, TAN, TDEN, S, C, XJ = 0;
    for (J = 0; J < NDIM * NDIM; J++) U[J] = ZER;
    for (I = 0; I < NDIM; I++) U[I * NDIM + I] = ONE;
iteration:
    ITCOUNT++;
    AMAX = ZER;
    for (I = 1; I < N; I++) {
        II = I * NDIM + I;
        JUP = I;
        for (J = 0; J < JUP; J++) {
            JJ = J * NDIM + J;
            IJ = J * NDIM + I;
            JI = I * NDIM + J;
            AII = A[II];
            AJJ = A[JJ];
            AOD = A[IJ];
            ASQ = AOD * AOD;
            if (ASQ > AMAX) AMAX = ASQ;
            if (ASQ <= EPS) continue;
            DIFFR = AII - AJJ;
            if (DIFFR < 0) {
                SIGN = -TWO;
                DIFFR = -DIFFR;
            } else {
                SIGN = TWO;
            }
            TDEN = DIFFR + sqrt(DIFFR * DIFFR + FOR * ASQ);
            TAN = SIGN * AOD / TDEN;
            C = ONE / sqrt(ONE + TAN * TAN);
            S = C * TAN;
            for (K = 0; K < N; K++) {
                KI = I * NDIM + K;
                KJ = J * NDIM + K;
                IK = K * NDIM + I;
                JK = K * NDIM + J;
                XJ = C * U[KJ] - S * U[KI];
                U[KI] = S * U[KJ] + C * U[KI];
                U[KJ] = XJ;
                if (K == J || K == I) continue;
                XJ = C * A[KJ] - S * A[KI];
                A[KI] = S * A[KJ] + C * A[KI];
                A[KJ] = XJ;
                A[IK] = A[KI];
                A[JK] = A[KJ];
            }
            A[II] = C * C * AII + S * S * AJJ + TWO * S * C * AOD;
            A[IJ] = ZER;
            A[JJ] = C * C * AJJ + S * S * AII - TWO * S * C * AOD;
            A[JI] = ZER;
        }
    }
    if (AMAX > EPS && ITCOUNT < MAXIT) goto iteration;
    for (I = 0; I < N; I++) {
        E[I] = A[I * NDIM + I];
    }
    if (AMAX > EPS) return 1;
    return 0;
}

int Jacobi_f(float* A, float* U, float* E, int NDIM, int N) {
    /*********************************************************************
     *   ARGUMENTS.                                                       *
     *        A       MATRIX TO BE DIAGONALIZED.                          *
     *        U       MATRIX OF EIGENVECTORS.                             *
     *        E       LIST OF EIGENVALUES                                 *
     *       NV       # OF EIGENVALUES NEEDED                             *
     *     NDIM       DIMENSION OF THE MATRICES AND THE VECTOR.           *
     *   Return       0: GOOD, 0: FAIL
     **********************************************************************/
    int MAXIT = 100, ITCOUNT = 0, JUP, I, J, K, II, JJ, KI, KJ, IJ, JI, IK, JK;
    float ZER = 0, ONE = 1, TWO = 2, FOR = 4, EPS = 1.0E-30;
    float AMAX, AII, AJJ, AOD, ASQ, DIFFR, SIGN, TAN, TDEN, S, C, XJ = 0;
    for (J = 0; J < NDIM * NDIM; J++) U[J] = ZER;
    for (I = 0; I < NDIM; I++) U[I * NDIM + I] = ONE;
iteration:
    ITCOUNT++;
    AMAX = ZER;
    for (I = 1; I < N; I++) {
        II = I * NDIM + I;
        JUP = I;
        for (J = 0; J < JUP; J++) {
            JJ = J * NDIM + J;
            IJ = J * NDIM + I;
            JI = I * NDIM + J;
            AII = A[II];
            AJJ = A[JJ];
            AOD = A[IJ];
            ASQ = AOD * AOD;
            if (ASQ > AMAX) AMAX = ASQ;
            if (ASQ <= EPS) continue;
            DIFFR = AII - AJJ;
            if (DIFFR < 0) {
                SIGN = -TWO;
                DIFFR = -DIFFR;
            } else {
                SIGN = TWO;
            }
            TDEN = DIFFR + sqrt(DIFFR * DIFFR + FOR * ASQ);
            TAN = SIGN * AOD / TDEN;
            C = ONE / sqrt(ONE + TAN * TAN);
            S = C * TAN;
            for (K = 0; K < N; K++) {
                KI = I * NDIM + K;
                KJ = J * NDIM + K;
                IK = K * NDIM + I;
                JK = K * NDIM + J;
                XJ = C * U[KJ] - S * U[KI];
                U[KI] = S * U[KJ] + C * U[KI];
                U[KJ] = XJ;
                if (K == J || K == I) continue;
                XJ = C * A[KJ] - S * A[KI];
                A[KI] = S * A[KJ] + C * A[KI];
                A[KJ] = XJ;
                A[IK] = A[KI];
                A[JK] = A[KJ];
            }
            A[II] = C * C * AII + S * S * AJJ + TWO * S * C * AOD;
            A[IJ] = ZER;
            A[JJ] = C * C * AJJ + S * S * AII - TWO * S * C * AOD;
            A[JI] = ZER;
        }
    }
    if (AMAX > EPS && ITCOUNT < MAXIT) goto iteration;
    for (I = 0; I < N; I++) {
        E[I] = A[I * NDIM + I];
    }
    if (AMAX > EPS) return 1;
    return 0;
}

void Eigen(double A[6][6], double U[6][6], double E[6]) {
    /*********************************************************************
     *   ARGUM6NTS.                                      *
     *        A       MATRIX TO BE DIAGONALIZED.                          *
     *        U       MATRIX OF EIGENVECTORS.                             *
     *        E       LIST OF EIGENVALUES                                 *
     *        N       # OF EIGENVALUES NEEDED                             *
     *        NDIM    DIMENSION OF THE MATRICES AND THE VECTOR.           *
     **********************************************************************/
    int NDIM = 6;
    int MAXIT = 100;
    int ITCOUNT = 0;
    double ZER = 0, ONE = 1, TWO = 2, FOR = 4, EPS = 1.0E-30;
    double AMAX, AII, AJJ, AOD, ASQ, DIFFR, SIGN, TAN, TDEN, S, C, XJ = 0;
    double AA[6][6];
    int JUP, I, J, K;
    for (J = 0; J < NDIM; J++) {
        for (I = 0; I < NDIM; I++) {
            U[I][J] = ZER;
            U[J][J] = ONE;
            AA[I][J] = A[I][J];
            AA[J][I] = A[J][I];
        }
    }
iteration:
    ITCOUNT++;
    AMAX = ZER;
    for (I = 1; I < NDIM; I++) {
        JUP = I;
        for (J = 0; J < JUP; J++) {
            AII = A[I][I];
            AJJ = A[J][J];
            AOD = A[I][J];
            ASQ = AOD * AOD;
            if (ASQ > AMAX) AMAX = ASQ;
            if (ASQ <= EPS) continue;
            DIFFR = AII - AJJ;
            if (DIFFR < 0) {
                SIGN = -TWO;
                DIFFR = -DIFFR;
            } else {
                SIGN = TWO;
            }
            TDEN = DIFFR + sqrt(DIFFR * DIFFR + FOR * ASQ);
            TAN = SIGN * AOD / TDEN;
            C = ONE / sqrt(ONE + TAN * TAN);
            S = C * TAN;
            for (K = 0; K < NDIM; K++) {
                XJ = C * U[K][J] - S * U[K][I];
                U[K][I] = S * U[K][J] + C * U[K][I];
                U[K][J] = XJ;
                if (K == J || K == I) continue;
                XJ = C * A[K][J] - S * A[K][I];
                A[K][I] = S * A[K][J] + C * A[K][I];
                A[K][J] = XJ;
                A[I][K] = A[K][I];
                A[J][K] = A[K][J];
            }
            A[I][I] = C * C * AII + S * S * AJJ + TWO * S * C * AOD;
            A[I][J] = ZER;
            A[J][J] = C * C * AJJ + S * S * AII - TWO * S * C * AOD;
            A[J][I] = ZER;
        }
    }
    if (AMAX > EPS && ITCOUNT < MAXIT) goto iteration;
    for (I = 0; I < NDIM; I++) {
        E[I] = A[I][I];
    }
    for (J = 0; J < NDIM; J++) {
        for (I = 0; I < NDIM; I++) {
            A[I][J] = AA[I][J];
            A[J][I] = AA[J][I];
        }
    }
}

/*

      KABSCH ALGORITHM FOR THE ALIGNMENT OF TWO RIGID STRUCTURES
      FIT X TO Y. Rotation matrix returned in U[3][3]. The return value
      of the call is the maximum distance

*/

double superpose(double* Y, double* X, int NATOMS, double U[][3]) {
    int I, J, K, I3;
    double R[3][3], S[3][3], RR[3][3], R1[9], VEC[9];
    double CX[3], CY[3], A[3][3], B[3][3], E[10];
    double EPS = 1.0E-15;
    double TWT, TEMP, TEMPX, TEMPY, TEMPZ, SUM, SCALE, PROD;
    double* WT = new double[NATOMS];
    double DIST, MAXDIST = 0;
    /*
     *     Reset origin to the center of sytems
     */
    for (I = 0; I < NATOMS; I++) {
        WT[I] = 1;
    }
    for (I = 0; I < 3; I++) {
        CX[I] = 0;
        CY[I] = 0;
    }
    TWT = 0;
    for (I = 0; I < NATOMS; I++) {
        I3 = I * 3;
        for (J = 0; J < 3; J++) {
            CX[J] = CX[J] + X[I3 + J] * WT[I];
            CY[J] = CY[J] + Y[I3 + J] * WT[I];
        }
        TWT = TWT + WT[I];
    }
    for (I = 0; I < 3; I++) {
        CX[I] = CX[I] / TWT;
        CY[I] = CY[I] / TWT;
    }
    for (I = 0; I < NATOMS; I++) {
        I3 = I * 3;
        for (J = 0; J < 3; J++) {
            X[I3 + J] = X[I3 + J] - CX[J];
            Y[I3 + J] = Y[I3 + J] - CY[J];
        }
    }
    /*
     *     UX->Y, i.e. Transfer X to fit Y
     */
    for (J = 0; J < 3; J++) {
        for (K = 0; K < 3; K++) {
            R[J][K] = 0;
            S[J][K] = 0;
        }
    }
    /*
     *     Build R and S matrices
     */
    for (I = 0; I < NATOMS; I++) {
        I3 = I * 3;
        for (J = 0; J < 3; J++) {
            for (K = 0; K < 3; K++) {
                R[J][K] = R[J][K] + WT[I] * Y[I3 + J] * X[I3 + K];
                S[J][K] = S[J][K] + WT[I] * X[I3 + J] * X[I3 + K];
            }
        }
    }
    /*
     *     Form RR matrix
     */
    for (J = 0; J < 3; J++) {
        for (K = 0; K < 3; K++) {
            RR[J][K] = 0;
            for (I = 0; I < 3; I++) {
                RR[J][K] = RR[J][K] + R[I][J] * R[I][K];
            }
        }
    }
    for (J = 0; J < 3; J++) {
        for (K = 0; K < 3; K++) {
            R1[K * 3 + J] = RR[J][K];
        }
    }
    Jacobi(R1, VEC, E, 3, 3);
    for (J = 0; J < 3; J++) {
        for (K = 0; K < 3; K++) {
            A[J][K] = VEC[K * 3 + J];
        }
    }

    /*
     *   Sort the three eigen values and the corresponding eigenvectors
     *
     *     i.e. E(0) <= E(1) <= E(2)
     */
    if (E[1] > E[0]) {
        TEMP = E[0];
        E[0] = E[1];
        E[1] = TEMP;
        for (I = 0; I < 3; I++) {
            TEMP = A[I][0];
            A[I][0] = A[I][1];
            A[I][1] = TEMP;
        }
    }
    if (E[2] > E[1]) {
        TEMP = E[1];
        E[1] = E[2];
        E[2] = TEMP;
        for (I = 0; I < 3; I++) {
            TEMP = A[I][1];
            A[I][1] = A[I][2];
            A[I][2] = TEMP;
        }
    }
    if (E[1] > E[0]) {
        TEMP = E[0];
        E[0] = E[1];
        E[1] = TEMP;
        for (I = 0; I < 3; I++) {
            TEMP = A[I][0];
            A[I][0] = A[I][1];
            A[I][1] = TEMP;
        }
    }
    /*
     *   CROSS PRODUCT: A3 = A1 x A2 (cross product), so that A1,A2,A3 form a
     *   right-handed system
     */
    A[0][2] = A[1][0] * A[2][1] - A[1][1] * A[2][0];
    A[1][2] = A[0][1] * A[2][0] - A[0][0] * A[2][1];
    A[2][2] = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    for (K = 0; K < 3; K++) {
        for (I = 0; I < 3; I++) {
            B[I][K] = 0;
            for (J = 0; J < 3; J++) {
                B[I][K] = B[I][K] + R[I][J] * A[J][K];
            }
        }
    }
    /*
     *   NORMALIZATION
     */
    if (E[0] < EPS) {
        B[0][0] = 1;
        B[1][0] = 0;
        B[2][0] = 0;
        B[0][1] = 0;
        B[1][1] = 1;
        B[2][1] = 0;
        B[0][2] = 0;
        B[1][2] = 0;
        B[2][2] = 1;
    } else if (E[1] < EPS) {
        /*
         *     ARBITARY VECTOR
         */
        B[0][1] = 1.0;
        B[1][1] = 3.14159265358979;
        B[2][1] = 2.0;
        SUM = 0.0;
        for (I = 0; I < 3; I++) {
            SUM = SUM + B[I][0] * B[I][0];
        }
        SCALE = 1 / sqrt(SUM);
        for (I = 0; I < 3; I++) {
            B[I][0] = B[I][0] * SCALE;
        }
        PROD = B[0][0] * B[0][1] + B[1][0] * B[1][1] + B[2][0] * B[2][1];
        B[0][1] = B[0][1] - PROD * B[0][0];
        B[1][1] = B[1][1] - PROD * B[1][0];
        B[2][1] = B[2][1] - PROD * B[2][0];
        SUM = 0;
        for (I = 0; I < 3; I++) {
            SUM = SUM + B[I][1] * B[I][1];
        }
        SCALE = 1 / sqrt(SUM);
        for (I = 0; I < 3; I++) {
            B[I][1] = B[I][1] * SCALE;
        }
    }
    for (K = 0; K < 2; K++) {
        SUM = 0;
        for (I = 0; I < 3; I++) {
            SUM = SUM + B[I][K] * B[I][K];
        }
        SCALE = 1 / sqrt(SUM);
        for (I = 0; I < 3; I++) {
            B[I][K] = B[I][K] * SCALE;
        }
    }
    /*
     *   CROSS PRODUCT
     */
    B[0][2] = B[1][0] * B[2][1] - B[1][1] * B[2][0];
    B[1][2] = B[0][1] * B[2][0] - B[0][0] * B[2][1];
    B[2][2] = B[0][0] * B[1][1] - B[0][1] * B[1][0];
    /*
     *   U = BA
     */
    for (I = 0; I < 3; I++) {
        for (J = 0; J < 3; J++) {
            U[I][J] = 0;
        }
    }
    for (K = 0; K < 3; K++) {
        for (I = 0; I < 3; I++) {
            for (J = 0; J < 3; J++) {
                U[I][J] = U[I][J] + B[I][K] * A[J][K];
            }
        }
    }
    /*
     *     X = UX
     */
    for (I = 0; I < NATOMS; I++) {
        TEMPX = 0;
        I3 = I * 3;
        for (J = 0; J < 3; J++) {
            TEMPX = TEMPX + U[0][J] * X[I3 + J];
        }
        TEMPY = 0;
        for (J = 0; J < 3; J++) {
            TEMPY = TEMPY + U[1][J] * X[I3 + J];
        }
        TEMPZ = 0;
        for (J = 0; J < 3; J++) {
            TEMPZ = TEMPZ + U[2][J] * X[I3 + J];
        }
        X[I3] = TEMPX;
        X[I3 + 1] = TEMPY;
        X[I3 + 2] = TEMPZ;
        DIST = sqrt((X[I3] - Y[I3]) * (X[I3] - Y[I3]) +
                    (X[I3 + 1] - Y[I3 + 1]) * (X[I3 + 1] - Y[I3 + 1]) +
                    (X[I3 + 2] - Y[I3 + 2]) * (X[I3 + 2] - Y[I3 + 2]));
        if (DIST > MAXDIST) {
            MAXDIST = DIST;
        }
    }
    delete[] WT;
    return MAXDIST;
}

//
// sort score array by descent, used in ROC plot
//
void sortScore(int sampSize, double* score, int* isGood) {
    int i, j;
    double tmpD;
    int tmpI;
    for (i = 0; i < sampSize; i++) {
        for (j = i + 1; j < sampSize; j++) {
            if (score[i] < score[j]) {
                tmpD = score[i];
                score[i] = score[j];
                score[j] = tmpD;

                tmpI = isGood[i];
                isGood[i] = isGood[j];
                isGood[j] = tmpI;
            }
        }
    }
}

//
// isGood: 1 is Good, 0 is not Good
// score is order by descent
//
// size of X and Y are 2+SampSize
void GenROCPoints(int SampSize, int* isGood, double* score, double* X, double* Y) {
    int i, j, k;
    int P = 0;
    int N = 0;
    int TP, FP;
    double TPrate, FPrate;

    for (i = 0; i < SampSize; i++) {
        if (isGood[i] == 1)
            P++;
        else
            N++;
    }

    X[0] = 0;
    Y[0] = 0;
    X[SampSize + 1] = 1;
    Y[SampSize + 1] = 1;

    TP = 0;
    FP = 0;
    for (i = 0; i < SampSize; i++) {
        if (1 == isGood[i])
            TP++;
        else
            FP++;

        TPrate = (double)TP / (double)P;
        FPrate = (double)FP / (double)N;
        X[i + 1] = FPrate;
        Y[i + 1] = TPrate;
    }
}

double CalAUC(int NPoints, double* X, double* Y) {
    double AUC = 0;
    double TrapzArea = 0;
    int i, j;

    j = NPoints - 1;
    for (i = 0; i < j; i++) {
        TrapzArea = (Y[i] + Y[i + 1]) * (X[i + 1] - X[i]) / 2;
        AUC += TrapzArea;
    }

    return AUC;
}

// generate enrichment factor
// opt: 0->1% , 1->5% , 2->10%, 3->20%
double calEF(int sampSize, int* isGood, double* score, int opt) {
    if (sampSize <= 0) {
        std::cerr << "no sample, sample size is 0\n";
        return 0;
    }

    int i;
    int P = 0;
    int n = 0;
    int TP;
    double TPrate, ActiveRate;

    if (0 == opt)
        n = sampSize * 0.01;
    else if (1 == opt)
        n = sampSize * 0.05;
    else if (2 == opt)
        n = sampSize * 0.1;
    else
        n = sampSize * 0.2;

    if (n <= 0) {
        std::cerr << "sample size is too small\n";
        return 0;
    }

    for (i = 0; i < sampSize; i++) {
        if (isGood[i] == 1) P++;
    }
    ActiveRate = (double)P / (double)sampSize;

    TP = 0;
    for (i = 0; i < n; i++) {
        if (1 == isGood[i]) TP++;
    }
    TPrate = (double)TP / (double)n;

    double EF = TPrate / ActiveRate;
    return EF;
}

// recursive call to genenrate combinations
void combinGen(int n, int m, int* iarray, int plevel, int* combines, int* np) {
    if (plevel == m) {
        for (int i = 0; i < m; i++) {
            combines[*np * m + i] = iarray[i];
        }
        (*np)++;
    } else {
        // Go to the next level
        int nplevel = plevel + 1;
        int* jarray = new int[n];
        //     for (int i=0; i < n; i++)
        //          jarray[i]=iarray[i];
        // combination
        for (int j = plevel; j < n; j++) {
            for (int i = 0; i < n; i++) jarray[i] = iarray[i];

            if (plevel > 0) {
                // before exchanging jarray[plevel] and jarray[j]
                // make sure jarray[j] > jarray[plevel-1]
                if (jarray[j] <= jarray[plevel - 1]) continue;
            }
            jarray[plevel] = iarray[j];
            jarray[j] = iarray[plevel];
            combinGen(n, m, jarray, nplevel, combines, np);
        }
        delete jarray;
    }
}
// Find the combination of m elements from n elements (n >= m)
int* combination(int n, int m) {
    int nx = 1;
    int mf = 1;
    int np = 0;
    int plevel = 0;
    for (int i = 0; i < m; i++) nx = nx * (n - i);
    for (int i = m; i > 0; i--) mf = mf * i;
    nx = nx / mf;

    int* combines = new int[nx * m];
    int* iarray = new int[n];
    for (int i = 0; i < n; i++) iarray[i] = i;
    combinGen(n, m, iarray, plevel, combines, &np);
    delete iarray;
    /*
       std::cout << "Final Permutations\n";
       std::cout << "number of combinations: " << np << "\n";
       for (int i=0; i < nx; i++) {
          for (int j=0; j < m; j++)
            std::cout << combines[i*m+j] << " ";
          std::cout << "\n";
       }
    */
    return combines;
}

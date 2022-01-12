#ifndef mathUtil_h
#define mathUtil_h
/******************************************************************************
 *                                                                             *
 *    Category:   Molecule Foundation Classes (MFC)                            *
 *    Function:   Implement math utilities for 3D operations, eigen value etc. *
 *    Author:     James Li                                                     *
 *    Date:       March, 2011                                                  *
 *                                                                             *
 ******************************************************************************/
#include <algorithm>
#include <string>
#include <vector>

void Rotation(double xyz0[3], double r[3][3], double xyz1[3]);
void RotMat(double alpha, double beta, double gama, double r[3][3]);

void Rotation_float(float xyz0[3], float r[3][3], float xyz1[3]);
void RotMat_float(float alpha, float beta, float gama, float r[3][3]);

void Euler2Mat(double alpha, double beta, double gama, double R[3][3], double DRA[3][3],
               double DRB[3][3], double DRG[3][3], double DDRAB[3][3], double DDRAG[3][3],
               double DDRBG[3][3], double DDRAA[3][3], double DDRBB[3][3], double DDRGG[3][3]);

// MathUtil
// Use a quaternion to represent a 3D rotation (axis and angle)
struct Quaternion {
    Quaternion();
    Quaternion(Quaternion* q);
    Quaternion(double angle, double* v);
    double x, y, z, w;
};

// Convert a quaternion (4 parameters) to a rotation matrix
void QuatToRMatrix(double m[3][3], double q[4]);

// Convert a quaternion object to a rotation matrix
void QuatToRMatrix(double m[3][3], Quaternion* Q);

// Rotation matrix to quaternion
void RMatrixToQuat(double m[3][3], Quaternion* Q);

// Quaternion multiplication (using 4-parameter form)
void QuatMultiple(double p1[4], double p2[4], double dest[4]);

// Quaternion multiplication (using quaternion objects)
void QuatMultiple(Quaternion* Q1, Quaternion* Q2, Quaternion* Q3);

// to judge if two Quaternions are equal
bool QuatEqual(Quaternion* Q1, Quaternion* Q2);

// Eigenvalue problem solver
void Eigen(double A[6][6], double U[6][6], double E[6]);

// Jacobi method for solving eigenvalue-eigenvector problems
int Jacobi(double* A, double* U, double* E, int NDIM, int NV);
int Jacobi_f(float* A, float* U, float* E, int NDIM, int NV);

// Kabsch algorithm for alignment of two rigid structure X to Y
// The rotation matrix also stored in rot[3][3] when return.
// The return value is the maximum distance of the matching points.
double superpose(double* Y, double* X, int NATOMS, double rot[][3]);

struct ROCRecord {
    int active;  // 1 is active, 0 is inactive
    double score;
};

// decendent order
struct compROCRecord : public std::binary_function<ROCRecord*, ROCRecord*, bool> {
    bool operator()(ROCRecord* A, ROCRecord* B) { return (A->score > B->score); }
};

struct scoreRecord {
    std::size_t idx;  // molecule index, start from 0
    double score;
    int active;
};

struct scoreRecord_SS {
    double corescore;
    double sidescore;
    int active;
    std::string name;
};

// decendent order
struct compScoreRecord : public std::binary_function<scoreRecord*, scoreRecord*, bool> {
    bool operator()(scoreRecord* A, scoreRecord* B) { return (A->score > B->score); }
};
//
// find max simi in confs of one database molecule
//
std::vector<scoreRecord> find_score_record(
    float* simi, const std::vector<uint16_t>& conformations_number_per_mol);

// sort score array by descent, isGood array changed with score
void sortScore(int sampSize, double* score, int* isGood);

// lixm.
// decendent order
struct compScoreRecord_SS : public std::binary_function<scoreRecord_SS*, scoreRecord_SS*, bool> {
    compScoreRecord_SS(double val = 0) : judge(val) {}
    bool operator()(scoreRecord_SS* A, scoreRecord_SS* B) {
        if (A->corescore > B->corescore) {
            if ((A->corescore - B->corescore) >= judge)
                return true;
            else {
                if (A->sidescore > B->sidescore) return true;
            }
        } else {
            if ((B->corescore - A->corescore) < judge) {
                if (A->sidescore > B->sidescore) return true;
            }
        }
        return false;
        // return (A->corescore > B->corescore);
    }

private:
    double judge;
};

// generate points of ROC curve
void GenROCPoints(int NPoints, int* isGood, double* score, double* X, double* Y);

// Calculate AUC value from points in ROC curve
double CalAUC(int NPoints, double* X, double* Y);

// generate enrichment factor
// opt: 0->1% , 1->5% , 2->10%, 3->20%
double calEF(int sampSize, int* isGood, double* score, int opt);

// combinations
int* combination(int n, int m);
void combinGen(int n, int m, int* iarray, int plevel, int* combines, int* np);

#endif /* mathUtil_h */

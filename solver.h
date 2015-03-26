#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "ofMain.h"
#include "Hemesh.h"
#include "DiscreteExteriorCalculus.h"

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::ConjugateGradient<SpMat> Cg;
typedef Eigen::VectorXf Vectorf;
typedef Eigen::VectorXd Vectord;

extern float maxEdgeLength;
extern float maxEdgeLengthSq;
extern vector<float> distances;

struct Force {
	float restLength;
	float restAngle;
	float length;
	float area;
	float invLength;
	float strength;
	ofVec3f dir;
	int index1,index2;
	int i0,i1,i2,i3;
	float bend, dBend, d2Bend;
};

void setupSolver(hemesh & mesh);
void geodesicDistance(hemesh & hmesh, vector<float> &out, bool tag=false);

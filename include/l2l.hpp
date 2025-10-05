#pragma once
#include <time.h>
#include <queue>
#include <vector>
#include "KDTree.hpp"
#include "line_fitting.hpp"
#include "line_3d.hpp"
#include "pmc.h"
#include <fstream>

#define PI 3.1415926536
#define SQRT3 1.732050808

typedef struct _ROTNODE
{
	float a, b, c, w;
	float ub, lb;
	int l;
	Eigen::Matrix3d SO3;
	friend bool operator < (const struct _ROTNODE & n1, const struct _ROTNODE & n2)
	{
		if(n1.ub != n2.ub)
			return n1.ub < n2.ub;
		else
			return n1.w < n2.w;
			//return n1.ub > n2.ub;
	}
	
}ROTNODE;

typedef struct _TRANSNODE
{
	float x, y, z, w;
	float ub, lb;
	int l;
	friend bool operator < (const struct _TRANSNODE & n1, const struct _TRANSNODE & n2)
	{
		if(n1.lb != n2.lb)
			return n1.lb > n2.lb;
		else
			return n1.w < n2.w;
			// return n1.ub > n2.ub;
	}
}TRANSNODE;

typedef struct _TRANSNODE_MC
{
	float x, y, z, w;
	float ub, lb;
	int l;
	friend bool operator < (const struct _TRANSNODE_MC & n1, const struct _TRANSNODE_MC & n2)
	{
		if(n1.ub != n2.ub)
			return n1.ub < n2.ub;
		else
			return n1.w < n2.w;
			// return n1.ub > n2.ub;
	}
}TRANSNODE_MC;


/********************************************************/



/********************************************************/

#define MAXROTLEVEL 5

class GoL2L
{
public:
	int Nm, Nd;
	std::vector<LINE3D> lModel, lData;
	std::vector<Eigen::Vector3d> d_Model, d_Data;
	std::vector<Eigen::Vector3d> b_Model, b_Data;

	ROTNODE initNodeRot;
	// RotNode initNodeRot;
	TRANSNODE initNodeTrans;
	TRANSNODE_MC initNodeTrans_MC;

	ROTNODE optNodeRot;
	// RotNode optNodeRot;
	TRANSNODE optNodeTrans;
	TRANSNODE_MC optNodeTrans_MC;

	optimizer optimizer_ = optimizer(10000, 0.000001);

	GoL2L();
	std::pair<double,int> Register();
	float RotationMSEThresh;
	float TranslationMSEThresh;
	float RotationSSEThresh;
	float TranslationSSEThresh;

	// float optError;
	int optInlier;
	Eigen::Matrix3d optR = Eigen::Matrix3d::Identity();
	Eigen::Vector3d optT = Eigen::Vector3d::Zero();

	clock_t clockBegin;
	int inlierNum;


private:
	//temp variables
	double * psi1;
	double * psi2;
	bool * is_positive_lb;
	std::vector<double> maxRotAng;
	Eigen::Matrix3d ToSkewSymmetric(const Eigen::Vector3d& rotation_vector);
	double ComputeTightRotBound(const ROTNODE & rotation_node, const Eigen::Vector3d& v,  const std::vector<Eigen::Vector3d> & vertex_rotation_vectors);
	std::vector<std::pair<int,int>> RotationSearch();
	void TranslationSearch(std::vector<std::pair<int,int>> Correspondence, std::vector<int> MaxClique);
	void TranslationMCSearch(std::vector<std::pair<int,int>> correspondence);
	std::vector<int> FindMaxClique(const std::vector<std::pair<int,int>> & correspondence);
	void Initialize();
	double ComputeDispDist(Eigen::Vector3d disp1, Eigen::Vector3d disp2);
	double ComputeLineDist(LINE3D l1, LINE3D l2);
	void Clear();
	
};

/********************************************************/

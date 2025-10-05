#pragma once
#include <time.h>
#include <queue>
#include <vector>
#include "p2p_fitting.hpp"
#include "plane_3d.hpp"
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

typedef struct _ROTNODE_FULL
{
	float a, b, c, w;
	float ub, lb;
	int l;
	Eigen::Matrix3d SO3;

	friend bool operator < (const struct _ROTNODE_FULL & n1, const struct _ROTNODE_FULL & n2)
	{
		if(n1.lb != n2.lb)
			return n1.lb > n2.lb;
		else
			return n1.w < n2.w;
			// return n1.ub > n2.ub;
	}
}ROTNODE_FULL;


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



/********************************************************/



/********************************************************/

#define MAXROTLEVEL 5

class GoP2P
{
public:
	int Nm, Nd;
	std::vector<PLANE3D> pModel, pData;
	std::vector<Eigen::Vector3d> d_Model, d_Data;
	std::vector<Eigen::Vector3d> n_Model, n_Data;

	ROTNODE initNodeRot;
	// RotNode initNodeRot;
	ROTNODE_FULL initNodeRot_full;

	TRANSNODE initNodeTrans;

	ROTNODE optNodeRot;
	ROTNODE_FULL optNodeRot_full;

	// RotNode optNodeRot;
	TRANSNODE optNodeTrans;
	plane_optimizer optimizer_ = plane_optimizer(10000, 0.000001);


	GoP2P();
	std::pair<double,int> Register();
	float RotationMSEThresh;
	float TranslationMSEThresh;
	float RotationSSEThresh;
	float TranslationSSEThresh;


	float optError;
	int optInlier;
	Eigen::Matrix3d optR = Eigen::Matrix3d::Identity();
	Eigen::Vector3d optT = Eigen::Vector3d::Zero();

	clock_t clockBegin;
	int inlierNum;


private:
	//temp variables
	double * psi1 = NULL;
	std::vector<double> maxRotAng;
	double * psi_rot_full;

	std::vector<int> RotationSearch();
	double RotationSearch_full();
	void TranslationSearch(std::vector<int> inlier_from_rot);
	void Initialize();
	double ComputeDispDist(Eigen::Vector3d disp1, Eigen::Vector3d disp2);
	void Clear();
	
	double OuterBnB();
	double InnerBnB(double * current_r_uncert, TRANSNODE* nodeTransOut, int level, Eigen::Matrix3d R0, double cur_maxRotAng);


};

/********************************************************/

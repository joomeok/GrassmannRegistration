#pragma once
#include <queue>
#include <vector>
#include <map>
#include <stdlib.h>
#include <math.h>
#include <chrono>
#include <vector>
#include <fstream>

#include "l2p_fitting.hpp"
#include "line_3d.hpp"
#include "plane_3d.hpp"


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

class GoL2P
{
public:
	int Nm, Nd;

    // Line Info
	std::map<int,std::vector<LINE3D>> lModel;
    // Plane Info
    std::vector<PLANE3D> pData;

	ROTNODE initNodeRot;
	ROTNODE_FULL initNodeRot_full;
	TRANSNODE initNodeTrans;
	TRANSNODE_MC initNodeTrans_MC;

	ROTNODE optNodeRot;
	ROTNODE_FULL optNodeRot_full;
	TRANSNODE optNodeTrans;
	TRANSNODE_MC optNodeTrans_MC;

	optimizer optimizer_ = optimizer(10000, 0.000001);


	GoL2P();
	std::pair<double,int> Register();
	float RotationMSEThresh;
	float TranslationMSEThresh;
	float RotationSSEThresh;
	float TranslationSSEThresh;

	// float optError;
	int optInlier_;
	Eigen::Matrix3d optR = Eigen::Matrix3d::Identity();
	Eigen::Vector3d optT = Eigen::Vector3d::Zero();

	clock_t clockBegin;
	int inlierNum;


private:
	//temp variables
	std::vector<double> * psi1;
	std::vector<double> maxRotAng;
	double * psi_rot_full;

	std::map<int,std::vector<LINE3D>> RotationSearch();
	double RotationSearch_full();
	double TranslationSearch(std::map<int,std::vector<LINE3D>> inlier_from_rot);
	void TranslationMCSearch(std::map<int,std::vector<LINE3D>> inlier_from_rot);

	void Initialize();

	void Clear();

};
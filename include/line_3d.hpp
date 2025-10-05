#pragma once
#include <eigen3/Eigen/Core>
#include "util_funcs.hpp"

class LINE3D{
	public:
	Eigen::Vector3d m_;
	Eigen::Vector3d d_;
	Eigen::Vector3d sp_;
	Eigen::Vector3d ep_;
	Eigen::Vector3d b_;
	bool empty_;
	int plane_id = 0;
	LINE3D();
	LINE3D(const Eigen::Vector3d & sp, const Eigen::Vector3d & ep);
	LINE3D(const Eigen::Vector3d & sp, const Eigen::Vector3d & ep, const int & id);

};

// double ComputeLineDistance(LINE3D l1, LINE3D l2);
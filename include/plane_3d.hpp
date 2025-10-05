#pragma once
#include <eigen3/Eigen/Core>
#include "util_funcs.hpp"

class PLANE3D{
	public:
	Eigen::Vector3d normal_;
	Eigen::Vector3d disp_;
    // Three points contained in the plane
	Eigen::Vector3d x1_;
	Eigen::Vector3d x2_;
    Eigen::Vector3d x3_;
    Eigen::Vector3d u_, v_;
    int id_;
	bool empty_;
	PLANE3D();
	PLANE3D(const Eigen::Vector3d & x1, const Eigen::Vector3d & x2, const Eigen::Vector3d & x3, const int & id);
};


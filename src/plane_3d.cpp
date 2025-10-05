#include "plane_3d.hpp"
#include <iostream>



PLANE3D::PLANE3D(){
	disp_ = Eigen::Vector3d::Zero();
	normal_ = Eigen::Vector3d::Zero();
}

PLANE3D::PLANE3D(const Eigen::Vector3d & x1, const Eigen::Vector3d & x2, const Eigen::Vector3d & x3, const int & id) : x1_(x1), x2_(x2), x3_(x3), id_(id) {
    normal_ = cross(x1_-x2_, x1_-x3_).normalized();
    // Ensure normal is directed as origin -> plane
    if(normal_.dot(x1_) < 0) normal_ = -normal_;
	// Set x1-x2 as the first basis
	u_ = (x1_-x2_).normalized();	
	v_ = cross(u_,normal_).normalized();

	disp_ = normal_.dot(x1_) * normal_;
}

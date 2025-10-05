#include "util_funcs.hpp"
#include <iostream>

double safe_acos(double value) {
    if (value<=-1.0) {
        return M_PI;
    } else if (value>=1.0) {
        return 0;
    } else {
        return std::acos(value);
    }
}

double ComputeGrassDist(Eigen::Vector4d v1, Eigen::Vector4d v2){	
	double acos;
    if(abs(v1.norm() - 1) > 1e-6 || abs(v2.norm() - 1) > 1e-6){
        std::cout << "Value Error: Norm is not 1" << std::endl;
        return -1;
    } 
	v1.dot(v2) > 0 ? acos = safe_acos(v1.dot(v2)) : acos = M_PI - safe_acos(v1.dot(v2));
	return acos;
}

Eigen::Vector3d cross(const Eigen::Vector3d & v1, const Eigen::Vector3d & v2){
	Eigen::Vector3d result;
	result << v1(1) * v2(2) - v1(2) * v2(1), v1(2) * v2(0) - v1(0) * v2(2), v1(0) * v2(1) - v1(1) * v2(0);
	return result;
}
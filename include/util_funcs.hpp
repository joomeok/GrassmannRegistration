#pragma once
#include <eigen3/Eigen/Core>

double safe_acos(double value);
double ComputeGrassDist(Eigen::Vector4d v1, Eigen::Vector4d v2);
Eigen::Vector3d cross(const Eigen::Vector3d & v1, const Eigen::Vector3d & v2);
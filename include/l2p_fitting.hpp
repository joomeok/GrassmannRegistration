#pragma once

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include "ceres/manifold.h"
#include "ceres/loss_function.h"
#include <sophus/se3.hpp>

#include "line_3d.hpp"
#include "plane_3d.hpp"

#include <cmath>
#include <eigen3/Eigen/Core>

class optimizer{
    public:
        optimizer(int max_iter_, double error_diff_);
    	int optimize_R(double RotationMSEThresh, std::map<int,std::vector<LINE3D>> lines, const std::vector<PLANE3D> & planes, const Eigen::Matrix3d & R_init,Eigen::Matrix3d & R_opt, int optInlier);
        double optimize_R_full(double RotationMSEThresh, std::map<int,std::vector<LINE3D>> lines, const std::vector<PLANE3D> & planes, const Eigen::Matrix3d & R_init,Eigen::Matrix3d & R_opt, double optErrorR, int Nm);

        double optimize_t(double TranslationMSEThresh, std::map<int,std::vector<LINE3D>>  lines, const std::vector<PLANE3D> & planes, std::map<int,std::vector<Eigen::Vector4d>> q2_fixed, const std::vector<int> & plane_id, const int & line_num, const Eigen::Vector3d & t_init, Eigen::Vector3d & t_opt, const Eigen::Matrix3d & R_opt, double optErrorT);
        int optimize_t_MC(double TranslationMSEThresh, std::map<int,std::vector<LINE3D>> lines, const std::vector<PLANE3D> & planes, std::map<int,std::vector<Eigen::Vector4d>> q2_fixed, 
const std::vector<int> & plane_id, const Eigen::Vector3d & t_init, Eigen::Vector3d & t_opt, const Eigen::Matrix3d & R_opt, int optInlier);
        int max_iter;
	    double error_diff;
};


class LineDirectionFunctor {
    public : 
        LineDirectionFunctor(const std::vector<double> & line_dir ,const std::vector<double> & u_vec,const std::vector<double> & v_vec);

    template <typename T>
    bool operator()(const T* const rotvec, T* residual) const;
    std::vector<double> line_dir_, u_, v_;
};

class LineDisplacementFunctor {
    public : 
        LineDisplacementFunctor(const std::vector<double> & line_disp, const std::vector<double> & fixed_q2, const std::vector<double> & plane_normal, const std::vector<double> & x, const Eigen::Matrix3d & optR);

    template <typename T>
    bool operator()(const T* const translation, T* residual) const;
    std::vector<double> b_,n_,x_, fixed_q2_;
    Eigen::Matrix3d optR_;
};
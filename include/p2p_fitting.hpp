#pragma once

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include "ceres/manifold.h"
#include "ceres/loss_function.h"
#include <sophus/se3.hpp>
#include "plane_3d.hpp"
#include <cmath>
#include <eigen3/Eigen/Core>

class plane_optimizer{
    public:
        plane_optimizer(int max_iter_, double error_diff_);

        // double optimize(const std::vector<PLANE3D> & data_planes, const Eigen::Matrix3d & R_init, const Eigen::Vector3d & t_init, Eigen::Matrix3d & R_opt, Eigen::Vector3d & t_opt);
    	// double optimize_R(const std::vector<PLANE3D> & data_planes, const Eigen::Matrix3d & R_init,Eigen::Matrix3d & R_opt, double MSEThreshold, int optInlier);
        double optimize_R_full(double RotationMSEThresh, const std::vector<PLANE3D> & data_planes, const std::vector<PLANE3D> & model_planes, const Eigen::Matrix3d & R_init,Eigen::Matrix3d & R_opt, double optErrorR, int Nm);
    	double optimize_R(const std::vector<PLANE3D> & data_planes,const std::vector<PLANE3D> & model_planes, const Eigen::Matrix3d & R_init,Eigen::Matrix3d & R_opt, double MSEThreshold, int optInlier);
        double optimize_t(const std::vector<PLANE3D> & data_planes, const std::vector<PLANE3D> & model_planes, const Eigen::Vector3d & t_init, Eigen::Vector3d & t_opt, const Eigen::Matrix3d & R_opt, double optErrorT);
        double optimize_R_t(const std::vector<PLANE3D> & data_planes, const std::vector<PLANE3D> & model_planes, const Eigen::Matrix3d & R_init,Eigen::Matrix3d & R_opt,const Eigen::Vector3d& t_init ,Eigen::Vector3d & t_opt, double optError);

        int max_iter;
	    double error_diff;
};

class PlaneDirectionFunctor {
    public : 
        PlaneDirectionFunctor(const std::vector<double> & model_p,const std::vector<double> & data_p);

    template <typename T>
    bool operator()(const T* const quaternion, T* residual) const;
    std::vector<double> model_p_, data_p_;
};

class PlaneDisplacementFunctor {
    public : 
        PlaneDisplacementFunctor(const std::vector<double> & model_l,const std::vector<double> & data_l,const Eigen::Matrix3d & optR);

    template <typename T>
    bool operator()(const T* const translation, T* residual) const;
    std::vector<double> model_p_disp_, data_p_;
    Eigen::Matrix3d optR_;
};



class PlaneOptimizeFunctor {
    public : 
        PlaneOptimizeFunctor(const std::vector<double> & data_n_x,const std::vector<double> & model_n_c);
    template <typename T>
    bool operator()(const T* const quaternion, const T* const translation, T* residual) const;
    std::vector<double> data_normal_x_, model_normal_c_;
};


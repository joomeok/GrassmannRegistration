#pragma once

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include "ceres/manifold.h"
#include "ceres/loss_function.h"
#include <sophus/se3.hpp>
#include "line_3d.hpp"
#include <cmath>
#include <eigen3/Eigen/Core>

class optimizer{
    public:
        optimizer(int max_iter_, double error_diff_);

        // double optimize(const std::vector<LINE3D> & data_lines, const Eigen::Matrix3d & R_init, const Eigen::Vector3d & t_init, Eigen::Matrix3d & R_opt, Eigen::Vector3d & t_opt);
    	// double optimize_R(const std::vector<LINE3D> & data_lines, const Eigen::Matrix3d & R_init,Eigen::Matrix3d & R_opt, double MSEThreshold, int optInlier);
    	double optimize_R(const std::vector<LINE3D> & data_lines,const std::vector<LINE3D> & model_lines, const Eigen::Matrix3d & R_init,Eigen::Matrix3d & R_opt, double MSEThreshold, int optInlier);
        double optimize_t(const std::vector<LINE3D> & data_lines, const std::vector<LINE3D> & model_lines, const Eigen::Vector3d & t_init, Eigen::Vector3d & t_opt, const Eigen::Matrix3d & R_opt);
        double optimize_t_MC(const std::vector<LINE3D> & data_lines, const std::vector<LINE3D> & model_lines, const Eigen::Vector3d & t_init, Eigen::Vector3d & t_opt, 
const Eigen::Matrix3d & R_opt, double MSEThreshold, int optInlier);
        int max_iter;
	    double error_diff;
};


// class FullPointLineFunctor {
//     public : 
//         FullPointLineFunctor(const std::vector<double> & model_l,const std::vector<double> & data_l);

//     template <typename T>
//     bool operator()(const T* const rotvec, const T* const translation, T* residual) const;
//     std::vector<double> model_l_, data_l_;
// };

class LineDirectionFunctor {
    public : 
        LineDirectionFunctor(const std::vector<double> & model_l,const std::vector<double> & data_l);

    template <typename T>
    bool operator()(const T* const rotvec, T* residual) const;
    std::vector<double> model_l_, data_l_;
};

class LineDisplacementFunctor {
    public : 
        LineDisplacementFunctor(const std::vector<double> & model_l,const std::vector<double> & data_l,const Eigen::Matrix3d & optR);

    template <typename T>
    bool operator()(const T* const translation, T* residual) const;
    std::vector<double> model_l_, data_l_;
    Eigen::Matrix3d optR_;
};

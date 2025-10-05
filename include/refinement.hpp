#pragma once 

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include "ceres/manifold.h"
#include "ceres/loss_function.h"
#include <cmath>
#include <eigen3/Eigen/Core>

class L2LFunctor {
    public : 
        L2LFunctor(const std::vector<double> & model_l,const std::vector<double> & data_l);

    template <typename T>
    bool operator()(const T* const rot, const T* const t,  T* residual) const;
    Eigen::Vector3d model_d_, data_d_, model_b_, data_b_;
};

class L2PFunctor {
    public : 
        L2PFunctor(const std::vector<double> & line, const std::vector<double> & plane);

    template <typename T>
    bool operator()(const T* const rot, const T* const t, T* residual) const;
    Eigen::Vector3d line_d_, line_b_, plane_u_, plane_v_, plane_c_;
};

class P2PFunctor {
    public : 
        P2PFunctor(const std::vector<double> & model_p,const std::vector<double> & data_p);

    template <typename T>
    bool operator()(const T* const rot, const T* const t, T* residual) const;
    Eigen::Vector3d model_u_, model_v_, data_u_, data_v_, model_c_, data_c_;
};

std::pair<Eigen::Matrix3d, Eigen::Vector3d> L2Loptimize(std::vector<std::vector<double>> lines1, const std::vector<std::vector<double>> & lines2, const Eigen::Matrix<double,4,4>& T_init);
std::pair<Eigen::Matrix3d, Eigen::Vector3d> L2Poptimize(std::vector<std::vector<double>> lines, const std::vector<std::vector<double>> & planes, const Eigen::Matrix<double,4,4> & T_init);
std::pair<Eigen::Matrix3d, Eigen::Vector3d> P2Poptimize(std::vector<std::vector<double>> planes1, const std::vector<std::vector<double>> & planes2, const Eigen::Matrix<double,4,4> & T_init);

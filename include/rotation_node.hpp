#pragma once
#include <eigen3/Eigen/Dense>

class RotNode
{   
    public:
        int n_side;
        int ring_idx;
        int phi_idx;
        int S1_idx;
        float ub, lb;
        int l;
        Eigen::Matrix3d SO3;
        bool is_Identity = false;
        friend bool operator < (const RotNode & n1, const RotNode & n2);
	
};
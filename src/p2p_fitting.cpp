#include "p2p_fitting.hpp"


PlaneDirectionFunctor::PlaneDirectionFunctor(const std::vector<double> & data_p, const std::vector<double> & model_p) : data_p_(data_p), model_p_(model_p){}

template <typename T>
bool PlaneDirectionFunctor::operator()(const T* const quaternion, T* residual) const {

    Eigen::Quaternion<T> q(quaternion[0], quaternion[1], quaternion[2], quaternion[3]);
    q.normalize();
    
    Eigen::Matrix<T, 3, 1> d_data;
    d_data << T(data_p_.at(0)), T(data_p_.at(1)), T(data_p_.at(2));
    Eigen::Matrix<T,3,1> d_R = q * d_data;
    Eigen::Matrix<T, 3, 1> d_model;
    d_model << T(model_p_.at(0)), T(model_p_.at(1)), T(model_p_.at(2));

    T dot1 = d_R(0,0) * T(model_p_.at(0)) + d_R(1,0) * T(model_p_.at(1)) + d_R(2,0) * T(model_p_.at(2));
    Eigen::Matrix<T,3,1> subtract;
    dot1 > 0 ? subtract = d_R- d_model : subtract = d_R + d_model;

    residual[0] = subtract[0];
    residual[1] = subtract[1];
    residual[2] = subtract[2];


    return true;
}

PlaneDisplacementFunctor::PlaneDisplacementFunctor(const std::vector<double> & data_p,const std::vector<double> & model_p_disp, const Eigen::Matrix3d & optR) : data_p_(data_p), model_p_disp_(model_p_disp), optR_(optR) {}

template <typename T>
bool PlaneDisplacementFunctor::operator()(const T* const translation, T* residual) const {

    Eigen::Matrix<T, 3, 1> t(translation[0], translation[1], translation[2]);
    Eigen::Matrix<T,3,1> model_disp;
    model_disp << T(model_p_disp_.at(0)), T(model_p_disp_.at(1)), T(model_p_disp_.at(2));
    
    Eigen::Matrix<T, 4, 1> model_disp_tilde;
    model_disp_tilde << T(model_p_disp_.at(0)), T(model_p_disp_.at(1)), T(model_p_disp_.at(2)) , T(1);
    model_disp_tilde.normalize();

    Eigen::Matrix<T, 3, 1> n_data;
    n_data << T(data_p_.at(6)), T(data_p_.at(7)), T(data_p_.at(8));


    Eigen::Matrix<T, 4, 1> Ru_data_bar;
    Ru_data_bar << T(data_p_.at(0)), T(data_p_.at(1)), T(data_p_.at(2)), T(0);

    Eigen::Matrix<T, 4, 1> Rv_data_bar;
    Rv_data_bar << T(data_p_.at(3)), T(data_p_.at(4)), T(data_p_.at(5)), T(0);

    Eigen::Matrix<T, 3, 1> x1_data;
    x1_data << T(data_p_.at(9)), T(data_p_.at(10)), T(data_p_.at(11));

    Eigen::Matrix<T, 3, 1> data_disp; 
    data_disp = n_data.dot(x1_data - t) * optR_ * n_data;

    Eigen::Matrix<T, 4, 1> data_disp_tilde;
    data_disp_tilde << data_disp(0,0), data_disp(1,0) , data_disp(2,0) , T(1);
    data_disp_tilde.normalize();

    T dot = data_disp_tilde(0,0) * model_disp_tilde(0,0) + data_disp_tilde(1,0) * model_disp_tilde(1,0) + data_disp_tilde(2,0) * model_disp_tilde(2,0) + data_disp_tilde(3,0) * model_disp_tilde(3,0);

    Eigen::Matrix<T, 4, 1> term1 = Ru_data_bar.dot(model_disp_tilde) * Ru_data_bar;
    Eigen::Matrix<T, 4, 1> term2 = Rv_data_bar.dot(model_disp_tilde) * Rv_data_bar;
    Eigen::Matrix<T, 4, 1> term3 = dot * data_disp_tilde;

    Eigen::Matrix<T, 4, 1> P_dm = term1 + term2 + term3;

    P_dm.normalize();



    Eigen::Matrix<T,4,1> subtract = P_dm - model_disp_tilde;
    residual[0] = subtract[0];
    residual[1] = subtract[1];
    residual[2] = subtract[2];
    residual[3] = subtract[3];

    return true;
}


PlaneOptimizeFunctor::PlaneOptimizeFunctor(const std::vector<double> & data_n_x,const std::vector<double> & model_n_c) : data_normal_x_(data_n_x), model_normal_c_(model_n_c) {}

template <typename T>
bool PlaneOptimizeFunctor::operator()(const T* const quaternion, const T* const translation, T* residual) const {


    Eigen::Quaternion<T> q(quaternion[0], quaternion[1], quaternion[2], quaternion[3]);
    q.normalize();
    Eigen::Matrix<T, 3, 1> d_data;
    d_data << T(data_normal_x_.at(0)), T(data_normal_x_.at(1)), T(data_normal_x_.at(2));
    Eigen::Matrix<T,3,1> d_R = q * d_data;

    T dot1 = d_R(0,0) * T(model_normal_c_.at(0)) + d_R(1,0) * T(model_normal_c_.at(1)) + d_R(2,0) * T(model_normal_c_.at(2));
    T acos1; 

    dot1 > 0 ? acos1 = acos(dot1) : acos1 = M_PI - acos(dot1);




    Eigen::Matrix<T, 3, 1> t(translation[0], translation[1], translation[2]);
    Eigen::Matrix<T,3,1> model_disp;
    model_disp << T(model_normal_c_.at(3)), T(model_normal_c_.at(4)), T(model_normal_c_.at(5));
    

    Eigen::Matrix<T, 3, 1> x1_data;
    x1_data << T(data_normal_x_.at(3)), T(data_normal_x_.at(4)), T(data_normal_x_.at(5));

    Eigen::Matrix<T, 3, 1> data_disp;
    T temp_dot = d_data.dot(x1_data -t); 
    data_disp = temp_dot *( q * d_data);

    T dot2 = model_disp.dot(data_disp);
    T acos2;
    dot2 > 0 ? acos2 = acos(dot2) : acos2 = M_PI - acos(dot2);

    residual[0] = sqrt(acos1 * acos1 + acos2 * acos2);

    return true;
}




plane_optimizer::plane_optimizer(int max_iter_, double error_diff_) : max_iter(max_iter_), error_diff(error_diff_){}

double plane_optimizer::optimize_R_t(const std::vector<PLANE3D> & data_planes, const std::vector<PLANE3D> & model_planes, const Eigen::Matrix3d & R_init,Eigen::Matrix3d & R_opt,const Eigen::Vector3d& t_init ,Eigen::Vector3d & t_opt, double optError){
    int plane_num = data_planes.size();
    double error = -1;
    Eigen::Matrix3d R_new = R_init;
    Eigen::Vector3d t_new = t_init;
    ceres::Problem problem;


    for(int iter = 0 ; iter < 5; iter++){
        ceres::Manifold* q_manifold = new ceres::EigenQuaternionManifold();
        Eigen::Quaterniond q_init(R_new);        
        Eigen::Matrix<double, 7, 1> initial_pose;
        initial_pose << q_init.w(), q_init.x(), q_init.y(), q_init.z(), t_new(0), t_new(1), t_new(2);
        problem.AddParameterBlock(initial_pose.data(), 4, q_manifold);  // For quaternion part
        problem.AddParameterBlock(initial_pose.data() + 4, 3);   

        for (int i = 0; i < plane_num; i++){
            Eigen::Vector3d d_R = R_new * data_planes.at(i).normal_;
            Eigen::Vector3d d_model = model_planes.at(i).normal_;
            Eigen::Vector3d c_R_t = (data_planes.at(i).normal_.dot(data_planes.at(i).x1_-t_new)) * R_new * data_planes.at(i).normal_;
            Eigen::Vector3d c_model = model_planes.at(i).disp_;

            double dot1 = d_R.dot(d_model);
            double temp_err1;
            dot1 > 0 ? temp_err1 = safe_acos(dot1) : temp_err1 = M_PI - safe_acos(dot1);
            double temp_err2 = ComputeGrassDist(c_R_t.homogeneous().normalized(), c_model.homogeneous().normalized());
            double temp_err = temp_err1 * temp_err1 + temp_err2 * temp_err2;

            std::vector<double> data_n_x = {data_planes.at(i).normal_(0), data_planes.at(i).normal_(1), data_planes.at(i).normal_(2), data_planes.at(i).x1_(0), data_planes.at(i).x1_(1), data_planes.at(i).x1_(2)};
            std::vector<double> model_n_c = {d_model(0), d_model(1), d_model(2), c_model(0),c_model(1), c_model(2)};  

            
            ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<PlaneOptimizeFunctor,1,4,3>(
            new PlaneOptimizeFunctor(data_n_x, model_n_c));  
            ceres::LossFunction* loss_function = new ceres::CauchyLoss(2.0); 
            problem.AddResidualBlock(cost_function, loss_function, initial_pose.data(), initial_pose.data() + 4);
        }

        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);

        Eigen::Quaterniond q_opt(initial_pose[0], initial_pose[1], initial_pose[2], initial_pose[3]);
        R_new = q_opt.normalized().toRotationMatrix();
        double updated_error = 0;
        t_new << initial_pose[4], initial_pose[5], initial_pose[6];


       for (int i = 0; i < plane_num; i++){
            Eigen::Vector3d d_R = R_new * data_planes.at(i).normal_;
            Eigen::Vector3d d_model = model_planes.at(i).normal_;
            Eigen::Vector3d c_R_t = (data_planes.at(i).normal_.dot(data_planes.at(i).x1_-t_new)) * R_new * data_planes.at(i).normal_;
            Eigen::Vector3d c_model = model_planes.at(i).disp_;

            double dot1 = d_R.dot(d_model);
            double temp_err1;
            dot1 > 0 ? temp_err1 = safe_acos(dot1) : temp_err1 = M_PI - safe_acos(dot1);
            double temp_err2 = ComputeGrassDist(c_R_t.homogeneous().normalized(), c_model.homogeneous().normalized());
            double temp_err = temp_err1 * temp_err1 + temp_err2 * temp_err2;
            updated_error += pow(temp_err,2);           

        }

       if(updated_error > 0 && updated_error - error < error_diff * plane_num){
			break;
        }
        error = updated_error;

    }

    if(error < optError){
        R_opt = R_new;
        t_opt = t_new;
    }
    else{
        R_opt = R_init;
        t_opt = t_init;
    }

    return error;

}




double plane_optimizer::optimize_R_full(double RotationMSEThresh, const std::vector<PLANE3D> & data_planes, const std::vector<PLANE3D> & model_planes, const Eigen::Matrix3d & R_init,Eigen::Matrix3d & R_opt, double optErrorR, int Nm){
    int plane_num = data_planes.size();
    double error = -1;
    Eigen::Matrix3d R_new = R_init;
    ceres::Problem problem;

    for(int iter = 0 ; iter < 5; iter++){
        ceres::Manifold* q_manifold = new ceres::EigenQuaternionManifold();
        Eigen::Quaterniond q_init(R_new);        
        Eigen::Matrix<double, 4, 1> initial_pose;
        initial_pose << q_init.w(), q_init.x(), q_init.y(), q_init.z();
        problem.AddParameterBlock(initial_pose.data(),4, q_manifold);
   

    for (int i = 0; i < plane_num; i++){
            Eigen::Vector3d d_R = R_new * data_planes.at(i).normal_;
            Eigen::Vector3d d_model = model_planes.at(i).normal_;

            double dot = d_R.dot(d_model);
            double temp_err;
            dot > 0 ? temp_err = safe_acos(dot) : temp_err = M_PI - safe_acos(dot);
    
            std::vector<double> model_normal = {d_model(0), d_model(1), d_model(2)};
    

            std::vector<double> data_normal = {data_planes.at(i).normal_(0), data_planes.at(i).normal_(1), data_planes.at(i).normal_(2)};
            ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<PlaneDirectionFunctor,3,4>(
            new PlaneDirectionFunctor(data_normal, model_normal));  

            ceres::LossFunction* loss_function = new ceres::CauchyLoss(2.0); 
            problem.AddResidualBlock(cost_function, loss_function, initial_pose.data());
    }


    
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);

        Eigen::Quaterniond q_opt(initial_pose[0], initial_pose[1], initial_pose[2], initial_pose[3]);
        R_new = q_opt.normalized().toRotationMatrix();
        double updated_error = 0;


    for (int i = 0; i < plane_num; i++){
            Eigen::Vector3d d_R = R_new * data_planes.at(i).normal_;
            Eigen::Vector3d d_model = model_planes.at(i).normal_;
            double dot = d_R.dot(d_model);
            double temp_err;
            dot > 0 ? temp_err = safe_acos(dot) : temp_err = M_PI - safe_acos(dot);

            updated_error += pow(temp_err,2);           
        }

        Eigen::Matrix3d R_i;


       if(updated_error > 0 && updated_error - error < error_diff * Nm){
			break;
        }

        error = updated_error;
    }


    if(error < optErrorR){
        R_opt = R_new;
    }
    else{
        R_opt = R_init;
    }

    return error;
}






double plane_optimizer::optimize_R(const std::vector<PLANE3D> & data_Planes, const std::vector<PLANE3D> & model_Planes, const Eigen::Matrix3d & R_init, Eigen::Matrix3d & R_opt, double MSEThreshold, int optInlier){
    int Plane_num = data_Planes.size();
    double error = -1;
    Eigen::Matrix3d R_new = R_init;
    ceres::Problem problem;
    Eigen::Quaterniond q_init(R_new);
    Eigen::Matrix<double, 4, 1> initial_pose;
    int inlier_new = 0;

    initial_pose << q_init.w(), q_init.x(), q_init.y(), q_init.z();

    ceres::Manifold* q_manifold = new ceres::EigenQuaternionManifold();
    problem.AddParameterBlock(initial_pose.data(),4, q_manifold);
    for (int i = 0; i < Plane_num; i++){
            Eigen::Vector3d d_R = R_new * data_Planes.at(i).normal_;
            Eigen::Vector3d d_model = model_Planes.at(i).normal_;

            double dot = d_R.dot(d_model);
            double ith_err;
            dot > 0 ? ith_err = safe_acos(dot) : ith_err = M_PI - safe_acos(dot);
            std::vector<double> model_normal = {d_model(0), d_model(1), d_model(2)};
    
            if(ith_err < MSEThreshold) inlier_new++;
            std::vector<double> data_normal = {data_Planes.at(i).normal_(0), data_Planes.at(i).normal_(1), data_Planes.at(i).normal_(2)};
            ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<PlaneDirectionFunctor,1,4>(
            new PlaneDirectionFunctor(data_normal, model_normal));  

            ceres::LossFunction* loss_function = new ceres::CauchyLoss(2.0); 
            problem.AddResidualBlock(cost_function, loss_function, initial_pose.data());
    }

    if(inlier_new > optInlier){
        R_opt = R_init;
        return inlier_new;
    } 
    
    else{
        ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = false;
    ceres::Solver::Summary summary;
    // Solve

    ceres::Solve(options, &problem, &summary);

    Eigen::Quaterniond q_opt(initial_pose[0], initial_pose[1], initial_pose[2], initial_pose[3]);
    R_new = q_opt.normalized().toRotationMatrix();

    R_opt = R_new;
    int current_Inlier = 0;
    for(int i = 0; i < Plane_num; i++)
    {	
        std::vector<double> data_direction;
        Eigen::Vector3d pData_n = R_opt * data_Planes.at(i).normal_;
        Eigen::Vector3d pModel = model_Planes.at(i).normal_;
        double dot = pData_n.dot(pModel);
        double error;
        dot > 0 ? error = safe_acos(dot) : error = M_PI - safe_acos(dot);

        if(error < MSEThreshold) current_Inlier++;
    }

    if(current_Inlier > optInlier) return current_Inlier;
    else return optInlier;
    }
}

double ComputeDispDist(Eigen::Vector3d disp1, Eigen::Vector3d disp2){
		Eigen::Vector4d b1, b2;
		b1 << disp1(0), disp1(1), disp1(2), 1;
		b2 << disp2(0), disp2(1), disp2(2), 1;
		b1 = b1 / b1.norm();
		b2 = b2 / b2.norm();
		double acos;
		b1.dot(b2) > 0 ? acos = std::acos(b1.dot(b2)) : acos = M_PI - std::acos(b1.dot(b2));
		double distc = acos;
		return distc;
}

double plane_optimizer::optimize_t(const std::vector<PLANE3D> & data_Planes, const std::vector<PLANE3D> & model_Planes, const Eigen::Vector3d & t_init, Eigen::Vector3d & t_opt, const Eigen::Matrix3d & R_opt, double optErrorT){
    int Plane_num = data_Planes.size();
    double error = -1;
    Eigen::Vector3d t_new = t_init;

    double error_before_opt = 0;
    for (int i = 0; i < Plane_num; i++){
            Eigen::Vector4d d_i_tilde = (data_Planes.at(i).normal_.dot(data_Planes.at(i).x1_ - t_new) * R_opt * data_Planes.at(i).normal_).homogeneous().normalized();
            Eigen::Vector4d d_m_tilde = model_Planes.at(i).disp_.homogeneous().normalized();
            
            Eigen::Vector3d Ru = R_opt * data_Planes.at(i).u_;
            Eigen::Vector3d Rv = R_opt * data_Planes.at(i).v_;
            Eigen::Vector4d Ru_bar, Rv_bar;
            Ru_bar << Ru(0), Ru(1), Ru(2), 0;
            Rv_bar << Rv(0), Rv(1), Rv(2), 0;

            Eigen::Vector4d P_dm = ((Ru_bar.dot(d_m_tilde)) * Ru_bar + (Rv_bar.dot(d_m_tilde)) * Rv_bar + (d_i_tilde.dot(d_m_tilde)) * d_i_tilde);
            P_dm.normalize();

            double acos;
            P_dm.dot(d_m_tilde) > 0 ? acos = safe_acos(P_dm.dot(d_m_tilde)) : acos = M_PI - safe_acos(P_dm.dot(d_m_tilde));


            error_before_opt += acos;
    }
    
    for(int iter = 0; iter < 5; iter++)
    {      
        ceres::Problem problem;
        Eigen::Matrix<double, 3, 1> initial_pose;

        initial_pose << t_new;

        problem.AddParameterBlock(initial_pose.data(), 3);
        for (int i = 0; i < Plane_num; i++){
                Eigen::Vector3d Ru = R_opt * data_Planes.at(i).u_;
				Eigen::Vector3d Rv = R_opt * data_Planes.at(i).v_;

                std::vector<double> data_p = {Ru(0), Ru(1), Ru(2), Rv(0), Rv(1), Rv(2), data_Planes.at(i).normal_(0), data_Planes.at(i).normal_(1), data_Planes.at(i).normal_(2), data_Planes.at(i).x1_(0), data_Planes.at(i).x1_(1), data_Planes.at(i).x1_(2)};
                std::vector<double> model_p_disp = {model_Planes.at(i).disp_(0), model_Planes.at(i).disp_(1), model_Planes.at(i).disp_(2)};                
                
                ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<PlaneDisplacementFunctor,4,3>(
                new PlaneDisplacementFunctor(data_p, model_p_disp, R_opt));  
                ceres::LossFunction* loss_function = new ceres::CauchyLoss(2.0); 
                problem.AddResidualBlock(cost_function, loss_function, initial_pose.data());
        }
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        // Solve
        ceres::Solve(options, &problem, &summary);
        t_new = initial_pose;
        double updated_error = 0;
        double euclidean_sum = 0;
        for (int i = 0; i < Plane_num; i++){
                Eigen::Vector4d d_i_tilde = (data_Planes.at(i).normal_.dot(data_Planes.at(i).x1_ - t_new) * R_opt * data_Planes.at(i).normal_).homogeneous().normalized();
                Eigen::Vector4d d_m_tilde = model_Planes.at(i).disp_.homogeneous().normalized();
                
                Eigen::Vector3d Ru = R_opt * data_Planes.at(i).u_;
				Eigen::Vector3d Rv = R_opt * data_Planes.at(i).v_;
				Eigen::Vector4d Ru_bar, Rv_bar;
				Ru_bar << Ru(0), Ru(1), Ru(2), 0;
				Rv_bar << Rv(0), Rv(1), Rv(2), 0;

				Eigen::Vector4d P_dm = ((Ru_bar.dot(d_m_tilde)) * Ru_bar + (Rv_bar.dot(d_m_tilde)) * Rv_bar + (d_i_tilde.dot(d_m_tilde)) * d_i_tilde);
				P_dm.normalize();


                double dist = (P_dm - d_m_tilde).norm() * (P_dm - d_m_tilde).norm(); 

                updated_error += dist;
        }

        
        std::cout << "After t: " << t_new(0) << " " << t_new(1) << " " << t_new(2) << std::endl;
        std::cout << "Before opt: " << error_before_opt <<" After opt: " << updated_error << std::endl;
    
       if(updated_error > 0 && updated_error - error < error_diff * Plane_num){
			break;
        }
        error = updated_error;
    }
    
    if(error < optErrorT){
        t_opt = t_new;
    }
    else{
        t_opt = t_init;
    } 

    return error;
}


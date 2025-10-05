#include "line_fitting.hpp"


LineDirectionFunctor::LineDirectionFunctor(const std::vector<double> & data_l,const std::vector<double> & model_l) : data_l_(data_l), model_l_(model_l){}

template <typename T>
bool LineDirectionFunctor::operator()(const T* const quaternion, T* residual) const {

    Eigen::Quaternion<T> q(quaternion[0], quaternion[1], quaternion[2], quaternion[3]);
    q.normalize();
    
    Eigen::Matrix<T, 3, 1> d_data;
    d_data << T(data_l_.at(0)), T(data_l_.at(1)), T(data_l_.at(2));
    Eigen::Matrix<T,3,1> d_R = q * d_data;

    T dot1 = d_R(0,0) * T(model_l_.at(0)) + d_R(1,0) * T(model_l_.at(1)) + d_R(2,0) * T(model_l_.at(2));
    T acos1; 

    T distc = T(1) - abs(dot1);
    residual[0] = distc;

    return true;
}

LineDisplacementFunctor::LineDisplacementFunctor(const std::vector<double> & data_l,const std::vector<double> & model_l, const Eigen::Matrix3d & optR) : data_l_(data_l), model_l_(model_l), optR_(optR) {}

template <typename T>
bool LineDisplacementFunctor::operator()(const T* const translation, T* residual) const {

    Eigen::Matrix<T, 3, 1> t(translation[0], translation[1], translation[2]);
    Eigen::Matrix<T,3,1> b_data;
    b_data << T(data_l_.at(0)), T(data_l_.at(1)), T(data_l_.at(2));
    
    Eigen::Matrix<T, 3, 1> d_data;
    d_data << T(data_l_.at(3)), T(data_l_.at(4)), T(data_l_.at(5));
    
    Eigen::Matrix<T,3,1> t_cross_d;
    t_cross_d << t(1,0) * d_data(2,0) - t(2,0) * d_data(1,0), t(2,0) * d_data(0,0) - t(0,0) * d_data(2,0), t(0,0) * d_data(1,0) - t(1,0) * d_data(0,0);
    Eigen::Matrix<T,3,1> d_cross_t_cross_d;
    d_cross_t_cross_d << d_data(1,0) * t_cross_d(2,0) - d_data(2,0) * t_cross_d(1,0), d_data(2,0) * t_cross_d(0,0) - d_data(0,0) * t_cross_d(2,0), 
    d_data(0,0) * t_cross_d(1,0) - d_data(1,0) * t_cross_d(0,0);




    Eigen::Matrix<T,3,1> data_disp;
    data_disp = optR_ * (b_data - d_cross_t_cross_d);

    Eigen::Matrix<T,3,1> model_disp;
    model_disp << T(model_l_.at(0)), T(model_l_.at(1)), T(model_l_.at(2));
    Eigen::Matrix<T,3,1> subtract = data_disp - model_disp;

    residual[0] = subtract[0];
    residual[1] = subtract[1];
    residual[2] = subtract[2];
    return true;
}



optimizer::optimizer(int max_iter_, double error_diff_) : max_iter(max_iter_), error_diff(error_diff_){}

    
double optimizer::optimize_R(const std::vector<LINE3D> & data_lines, const std::vector<LINE3D> & model_lines, const Eigen::Matrix3d & R_init, Eigen::Matrix3d & R_opt, double MSEThreshold, int optInlier){
    int line_num = data_lines.size();
    double error = -1;
    Eigen::Matrix3d R_new = R_init;

    ceres::Problem problem;
    Eigen::Quaterniond q_init(R_new);
    Eigen::Matrix<double, 4, 1> initial_pose;
    int inlier_new = 0;

    initial_pose << q_init.w(), q_init.x(), q_init.y(), q_init.z();

    ceres::Manifold* q_manifold = new ceres::EigenQuaternionManifold();
    problem.AddParameterBlock(initial_pose.data(),4, q_manifold);
    for (int i = 0; i < line_num; i++){
            Eigen::Vector3d d_R = R_new * data_lines.at(i).d_;
            // std::vector<double> data_R = { d_R(0), d_R(1), d_R(2)};
            Eigen::Vector3d d_model = model_lines.at(i).d_;

            double dot = d_R.dot(d_model);
            double ith_err;
            dot > 0 ? ith_err = acos(dot) : ith_err = M_PI - acos(dot);


       
            std::vector<double> nearest_dir = {d_model(0), d_model(1), d_model(2)};


            if(ith_err < 0.008) inlier_new++;
            std::vector<double> data_dir = {data_lines.at(i).d_(0), data_lines.at(i).d_(1), data_lines.at(i).d_(2)};
            ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<LineDirectionFunctor,1,4>(
            new LineDirectionFunctor(data_dir, nearest_dir));  

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
    // std::cout << "R_opt : "  << R_opt << std::endl;
    int current_Inlier = 0;
    for(int i = 0; i < line_num; i++)
    {	
        std::vector<double> data_direction;
        Eigen::Vector3d lData_R = R_opt * data_lines.at(i).d_;
        Eigen::Vector3d lModel = model_lines.at(i).d_;

        double dot = lData_R.dot(lModel);
        double error;
        dot > 0 ? error = acos(dot) : error = M_PI - acos(dot);

        if(error < 0.008) current_Inlier++;
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
double optimizer::optimize_t(const std::vector<LINE3D> & data_lines, const std::vector<LINE3D> & model_lines, const Eigen::Vector3d & t_init, Eigen::Vector3d & t_opt, const Eigen::Matrix3d & R_opt){
    int line_num = data_lines.size();
    double error = -1;
    Eigen::Vector3d t_new = t_init;

    for(int iter = 0; iter < max_iter; iter++)
    {      
        ceres::Problem problem;
        Eigen::Matrix<double, 3, 1> initial_pose;

        initial_pose << t_new;
    

        problem.AddParameterBlock(initial_pose.data(), 3);
        // double tot_dist = 0;
        for (int i = 0; i < line_num; i++){
                std::vector<double> data_l = {data_lines.at(i).b_(0), data_lines.at(i).b_(1), data_lines.at(i).b_(2), data_lines.at(i).d_(0), data_lines.at(i).d_(1), data_lines.at(i).d_(2)};
                std::vector<double> model_l = {model_lines.at(i).b_(0), model_lines.at(i).b_(1), model_lines.at(i).b_(2), model_lines.at(i).d_(0), model_lines.at(i).d_(1), model_lines.at(i).d_(2)};                
                
                ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<LineDisplacementFunctor,3,3>(
                new LineDisplacementFunctor(data_l, model_l, R_opt));  

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
        for (int i = 0; i < line_num; i++){
                Eigen::Vector3d data_disp_Rt = R_opt * (data_lines.at(i).b_ - data_lines.at(i).d_.cross(t_new.cross(data_lines.at(i).d_)));
                Eigen::Vector3d model_disp = model_lines.at(i).b_;
                updated_error += pow(std::acos(data_disp_Rt.homogeneous().normalized().dot(model_disp.homogeneous().normalized())),2);
        }
        updated_error = sqrt(updated_error);

        
 
       if(updated_error > 0 && updated_error - error < error_diff * line_num){
			break;
        }
        error = updated_error;
    }


    t_opt = t_new;

    return error;
}


double optimizer::optimize_t_MC(const std::vector<LINE3D> & data_lines, const std::vector<LINE3D> & model_lines, const Eigen::Vector3d & t_init, Eigen::Vector3d & t_opt, 
const Eigen::Matrix3d & R_opt, double MSEThreshold, int optInlier){
    int line_num = data_lines.size();
    double error = -1;
    Eigen::Vector3d t_new = t_init;
 
    ceres::Problem problem;
    Eigen::Matrix<double, 3, 1> initial_pose;
    int inlier_new = 0;
    initial_pose << t_new;
    problem.AddParameterBlock(initial_pose.data(),3);
    std::vector<Eigen::Vector3d> matched_data_l;
    for (int i = 0; i < line_num; i++){
        Eigen::Vector3d data_disp_Rt = R_opt * (data_lines.at(i).b_ - data_lines.at(i).d_.cross(t_new.cross(data_lines.at(i).d_)));
        std::vector<double> data_disp_Rt_vec = {data_disp_Rt(0),data_disp_Rt(1),data_disp_Rt(2)};
        Eigen::Vector3d model_disp = model_lines.at(i).b_;
        
        double dot = data_disp_Rt.homogeneous().normalized().dot(model_disp.homogeneous().normalized());

        double ith_err;
        dot > 0 ? ith_err = acos(dot) : ith_err = M_PI - acos(dot);

      
        std::vector<double> data_l = {data_lines.at(i).b_(0), data_lines.at(i).b_(1), data_lines.at(i).b_(2), data_lines.at(i).d_(0), data_lines.at(i).d_(1), data_lines.at(i).d_(2)};
        std::vector<double> model_l = {model_lines.at(i).b_(0), model_lines.at(i).b_(1), model_lines.at(i).b_(2), model_lines.at(i).d_(0), model_lines.at(i).d_(1), model_lines.at(i).d_(2)};        
        if(ith_err < 0.008) inlier_new++;
        ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<LineDisplacementFunctor,3,3>(
        new LineDisplacementFunctor(data_l, model_l, R_opt));  
        ceres::LossFunction* loss_function = new ceres::CauchyLoss(2.0); 
        problem.AddResidualBlock(cost_function, loss_function, initial_pose.data());
        }

        if(inlier_new > optInlier){
            t_opt = t_init;
            return inlier_new;
        } 
        
        else{
            ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        // Solve

        ceres::Solve(options, &problem, &summary);

        t_opt = initial_pose;

    int current_Inlier = 0;

    for(int i = 0; i < line_num; i++)
    {	
        Eigen::Vector3d data_disp_Rt = R_opt * (data_lines.at(i).b_ - data_lines.at(i).d_.cross(t_opt.cross(data_lines.at(i).d_)));
        Eigen::Vector3d model_disp = model_lines.at(i).b_;

        double error = std::acos(data_disp_Rt.homogeneous().normalized().dot(model_disp.homogeneous().normalized()));


        if(error < 0.008) current_Inlier++;
    }
    if(current_Inlier > optInlier) return current_Inlier;
    else return optInlier;
    }
}


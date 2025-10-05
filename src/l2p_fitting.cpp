#include "l2p_fitting.hpp"

LineDirectionFunctor::LineDirectionFunctor(const std::vector<double> & line_dir ,const std::vector<double> & u_vec,const std::vector<double> & v_vec) : line_dir_(line_dir), u_(u_vec), v_(v_vec){}

template <typename T>
bool LineDirectionFunctor::operator()(const T* const quaternion, T* residual) const {

    Eigen::Quaternion<T> q(quaternion[0], quaternion[1], quaternion[2], quaternion[3]);
    q.normalize();
    
    Eigen::Matrix<T, 3, 1> d_, u_vec, v_vec, u_R, v_R;
    d_ << T(line_dir_.at(0)), T(line_dir_.at(1)), T(line_dir_.at(2));
    u_vec << T(u_.at(0)), T(u_.at(1)), T(u_.at(2));
    v_vec << T(v_.at(0)), T(v_.at(1)), T(v_.at(2));
    u_R = q * u_vec;
    v_R = q * v_vec;
    
    T d_dot_u_R = d_(0) * u_R(0) + d_(1) * u_R(1) + d_(2) * u_R(2);
    T d_dot_v_R = d_(0) * v_R(0) + d_(1) * v_R(1) + d_(2) * v_R(2);

    Eigen::Matrix<T,3,1> q_R = d_dot_u_R * u_R + d_dot_v_R * v_R;
    q_R = q_R / q_R.norm();
    T dot1 = q_R(0,0) * T(line_dir_.at(0)) + q_R(1,0) * T(line_dir_.at(1)) + q_R(2,0) * T(line_dir_.at(2));
    T acos1; 

    T distc = T(1) - abs(dot1);
    residual[0] = distc;
    return true;
}

LineDisplacementFunctor::LineDisplacementFunctor(const std::vector<double> & line_disp, const std::vector<double> & fixed_q2, const std::vector<double> & plane_normal, const std::vector<double> & x, const Eigen::Matrix3d & optR) : b_(line_disp),fixed_q2_(fixed_q2) ,n_(plane_normal), x_(x), optR_(optR) {}

template <typename T>
bool LineDisplacementFunctor::operator()(const T* const translation, T* residual) const {

    Eigen::Matrix<T, 3, 1> t(translation[0], translation[1], translation[2]);

    Eigen::Matrix<T,4,1> fixed_q2;
    fixed_q2 << T(fixed_q2_.at(0)), T(fixed_q2_.at(1)), T(fixed_q2_.at(2)), T(fixed_q2_.at(3));

    Eigen::Matrix<T,4,1> line_b;
    line_b << T(b_.at(0)), T(b_.at(1)), T(b_.at(2)), T(b_.at(3));
    
     
    Eigen::Matrix<T, 3, 1> plane_x;
    plane_x << T(x_.at(0)), T(x_.at(1)), T(x_.at(2));
    Eigen::Matrix<T, 3, 1> plane_n;
    plane_n << T(n_.at(0)), T(n_.at(1)), T(n_.at(2));

    Eigen::Matrix<T,3,1> c;

    c = (plane_n[0] * (plane_x - t)[0] + plane_n[1] * (plane_x - t)[1] + plane_n[2] * (plane_x - t)[2])* optR_* plane_n;
    
    Eigen::Matrix<T,4,1> c_tilde;


    c_tilde << c(0,0), c(1,0), c(2,0), T(1);
    c_tilde.normalize();

    Eigen::Matrix<T,4,1> q2 = fixed_q2 + line_b.dot(c_tilde) * c_tilde;
    // q2 = q2 / q2.norm();
    q2.normalize();

    Eigen::Matrix<T,4,1> subtract = line_b - q2;

    residual[0] = subtract[0];
    residual[1] = subtract[1];
    residual[2] = subtract[2];
    residual[3] = subtract[3];



    return true;
}


optimizer::optimizer(int max_iter_, double error_diff_) : max_iter(max_iter_), error_diff(error_diff_){} 

double optimizer::optimize_R_full(double RotationMSEThresh, std::map<int,std::vector<LINE3D>> lines, const std::vector<PLANE3D> & planes, const Eigen::Matrix3d & R_init,Eigen::Matrix3d & R_opt, double optErrorR, int Nm){
    int plane_num = planes.size();
    double error = -1;
    Eigen::Matrix3d R_new = R_init;
    ceres::Problem problem;



    for(int iter = 0 ; iter < max_iter; iter++){
        ceres::Manifold* q_manifold = new ceres::EigenQuaternionManifold();
        
        Eigen::Quaterniond q_init(R_new);        
        Eigen::Matrix<double, 4, 1> initial_pose;
        initial_pose << q_init.w(), q_init.x(), q_init.y(), q_init.z();
        problem.AddParameterBlock(initial_pose.data(),4, q_manifold);
   
        for (int i = 0; i < plane_num; i++){
            int line_num = lines[i].size();
            Eigen::Vector3d u = planes.at(i).u_;
            Eigen::Vector3d v = planes.at(i).v_;
            for(int j=0; j< line_num; j++){
                Eigen::Vector3d d = lines[i].at(j).d_;
                Eigen::Vector3d q1_R = ((d.dot(R_new * u)) * R_new * u + (d.dot(R_new * v)) * R_new * v).normalized();
                double dot = d.dot(q1_R);
                double temp_err;
                dot > 0 ? temp_err = acos(dot) : temp_err = M_PI - acos(dot);

                std::vector<double> line_dir = {d(0), d(1), d(2)};
                std::vector<double> u_vec = {u(0), u(1), u(2)};
                std::vector<double> v_vec = {v(0), v(1), v(2)};

                ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<LineDirectionFunctor,1,4>(
                new LineDirectionFunctor(line_dir, u_vec,v_vec));  
                ceres::LossFunction* loss_function = new ceres::CauchyLoss(2.0); 
                problem.AddResidualBlock(cost_function, loss_function, initial_pose.data());
            }
        }
    
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);

        Eigen::Quaterniond q_opt(initial_pose[0], initial_pose[1], initial_pose[2], initial_pose[3]);
        R_new = q_opt.normalized().toRotationMatrix();
        double updated_error = 0;


        for(int i = 0; i < plane_num; i++)
        {	
            int line_num = lines[i].size();
            Eigen::Vector3d u = planes.at(i).u_;
            Eigen::Vector3d v = planes.at(i).v_;

            for(int j = 0; j<line_num; j++){
                Eigen::Vector3d d = lines[i].at(j).d_;
                Eigen::Vector3d q1_R = ((d.dot(R_new * u)) * R_new * u + (d.dot(R_new * v)) * R_new * v).normalized();
                double dot = d.dot(q1_R);
                double temp_err;
                dot > 0 ? temp_err = acos(dot) : temp_err = M_PI - acos(dot);
                updated_error += pow(temp_err,2);
            }
        }


       if(updated_error > 0 && updated_error - error < error_diff * Nm){
			break;
        }
        std::cout << "Error after opt: " << updated_error << std::endl;

        error = updated_error;
    }


    if(error < optErrorR){

        R_opt = R_new;
        
    }
    else{
        R_opt = R_init;
    }
    std::cout << "After opt: " << R_opt << std::endl;

    return error;
}



int optimizer::optimize_R(double RotationMSEThresh, std::map<int,std::vector<LINE3D>> lines, const std::vector<PLANE3D> & planes, const Eigen::Matrix3d & R_init,Eigen::Matrix3d & R_opt, int optInlier){
    int plane_num = planes.size();
    double error = -1;
    Eigen::Matrix3d R_new = R_init;

    ceres::Problem problem;
    Eigen::Quaterniond q_init(R_new);
    Eigen::Matrix<double, 4, 1> initial_pose;
    int inlier_new = 0;
    
    initial_pose << q_init.w(), q_init.x(), q_init.y(), q_init.z();
    ceres::Manifold* q_manifold = new ceres::EigenQuaternionManifold();
    problem.AddParameterBlock(initial_pose.data(),4, q_manifold);
    for (int i = 0; i < plane_num; i++){
        int line_num = lines[i].size();
        Eigen::Vector3d u = planes.at(i).u_;
        Eigen::Vector3d v = planes.at(i).v_;
        for(int j=0; j< line_num; j++){
            Eigen::Vector3d d = lines[i].at(j).d_;
            Eigen::Vector3d q1_R = ((d.dot(R_new * u)) * R_new * u + (d.dot(R_new * v)) * R_new * v).normalized();
            double dot = d.dot(q1_R);
            double temp_err;
            dot > 0 ? temp_err = acos(dot) : temp_err = M_PI - acos(dot);
            if(temp_err < RotationMSEThresh) inlier_new++;
            std::vector<double> line_dir = {d(0), d(1), d(2)};
            std::vector<double> u_vec = {u(0), u(1), u(2)};
            std::vector<double> v_vec = {v(0), v(1), v(2)};

            ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<LineDirectionFunctor,1,4>(
            new LineDirectionFunctor(line_dir, u_vec,v_vec));  
            ceres::LossFunction* loss_function = new ceres::CauchyLoss(2.0); 
            problem.AddResidualBlock(cost_function, loss_function, initial_pose.data());

        }


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
    for(int i = 0; i < plane_num; i++)
    {	
        int line_num = lines[i].size();
        Eigen::Vector3d u = planes.at(i).u_;
        Eigen::Vector3d v = planes.at(i).v_;

        for(int j = 0; j<line_num; j++){
            Eigen::Vector3d d = lines[i].at(j).d_;
            Eigen::Vector3d q1_R = ((d.dot(R_new * u)) * R_new * u + (d.dot(R_new * v)) * R_new * v).normalized();
            double dot = d.dot(q1_R);
            double temp_err;
            dot > 0 ? temp_err = acos(dot) : temp_err = M_PI - acos(dot);
            if(temp_err < RotationMSEThresh) current_Inlier++;

        }
    }


    if(current_Inlier > optInlier) return current_Inlier;
    else return optInlier;
    }
}


double optimizer::optimize_t(double TranslationMSEThresh, std::map<int,std::vector<LINE3D>> lines, const std::vector<PLANE3D> & planes, std::map<int,std::vector<Eigen::Vector4d>> q2_fixed, 
const std::vector<int> & plane_id, const int & line_num, const Eigen::Vector3d & t_init, Eigen::Vector3d & t_opt, const Eigen::Matrix3d & R_opt, double optErrorT){
    int plane_num = plane_id.size();

    double error = -1;
    Eigen::Vector3d t_new = t_init;
    for(int iter = 0; iter < max_iter; iter++)
    {      
        ceres::Problem problem;
        Eigen::Matrix<double, 3, 1> initial_pose;

        initial_pose << t_new;
    
        problem.AddParameterBlock(initial_pose.data(), 3);
        // double tot_dist = 0;
        for (int i = 0; i < plane_num; i++){
                int id = plane_id.at(i);
                int line_num = lines[id].size();
                std::vector<double> n = {planes.at(id).normal_(0), planes.at(id).normal_(1), planes.at(id).normal_(2)};
                std::vector<double> x = {planes.at(id).x1_(0), planes.at(id).x1_(1), planes.at(id).x1_(2)};

                for (int j=0; j< line_num; j++){
                    Eigen::Vector4d b_tilde = lines[id].at(j).b_.homogeneous().normalized();      
                    std::vector<double> b_tilde_vec = {b_tilde(0),b_tilde(1), b_tilde(2), b_tilde(3)};

                    std::vector<double> curr_q2_fixed = {q2_fixed[id].at(j)(0), q2_fixed[id].at(j)(1), q2_fixed[id].at(j)(2), q2_fixed[id].at(j)(3)};

                    ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<LineDisplacementFunctor,4,3>(
                    new LineDisplacementFunctor(b_tilde_vec, curr_q2_fixed, n,x, R_opt));  

                    ceres::LossFunction* loss_function = new ceres::CauchyLoss(2.0); 
                    problem.AddResidualBlock(cost_function, loss_function, initial_pose.data());
                }
        }
    
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        // Solve
        ceres::Solve(options, &problem, &summary);
        t_new = initial_pose;

        double updated_error = 0;
        int current_Inlier = 0;

        for (int i = 0; i < plane_num; i++){
                int id = plane_id.at(i);
                int line_num = lines[id].size();
                Eigen::Vector3d n = planes.at(id).normal_;
                Eigen::Vector3d x = planes.at(id).x1_;
                Eigen::Vector4d c_new_tilde= (n.dot(x-t_new) * R_opt * n).homogeneous().normalized();

                for (int j=0; j< line_num; j++){
                    Eigen::Vector4d b_tilde = lines[id].at(j).b_.homogeneous().normalized();      

                    Eigen::Vector4d q2_new = (q2_fixed[id].at(j) + b_tilde.dot(c_new_tilde) * c_new_tilde).normalized();

                    double dot = b_tilde.dot(q2_new);
                    double temp_err;
                    dot > 0 ? temp_err = acos(dot) : temp_err = M_PI - acos(dot);
                    if(temp_err < TranslationMSEThresh) current_Inlier++;
                    updated_error += pow(ComputeGrassDist(b_tilde,q2_new),2);
                }
        }

        // updated_error = sqrt(updated_error);

        
       if(updated_error > 0 && updated_error - error < error_diff * line_num){
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


int optimizer::optimize_t_MC(double TranslationMSEThresh, std::map<int,std::vector<LINE3D>> lines, const std::vector<PLANE3D> & planes, std::map<int,std::vector<Eigen::Vector4d>> q2_fixed, 
const std::vector<int> & plane_id, const Eigen::Vector3d & t_init, Eigen::Vector3d & t_opt, const Eigen::Matrix3d & R_opt, int optInlier){
    int plane_num = plane_id.size();
    ceres::Problem problem;
    Eigen::Vector3d t_new = t_init;
    Eigen::Matrix<double, 3, 1> initial_pose;
    initial_pose << t_new;
    int inlier_new = -1;

    problem.AddParameterBlock(initial_pose.data(), 3);
    for (int i = 0; i < plane_num; i++){
            int id = plane_id.at(i);
            int line_num = lines[id].size();
            std::vector<double> n_vec = {planes.at(id).normal_(0), planes.at(id).normal_(1), planes.at(id).normal_(2)};
            std::vector<double> x_vec = {planes.at(id).x1_(0), planes.at(id).x1_(1), planes.at(id).x1_(2)};
            Eigen::Vector4d c_new_tilde= (planes.at(id).normal_.dot(planes.at(id).x1_-t_new) * R_opt * planes.at(id).normal_).homogeneous().normalized();

            for (int j=0; j< line_num; j++){
                Eigen::Vector4d b_tilde = lines[id].at(j).b_.homogeneous().normalized();      
                Eigen::Vector4d q2_new = (q2_fixed[id].at(j) + b_tilde.dot(c_new_tilde) * c_new_tilde).normalized();

                double dot = b_tilde.dot(q2_new);
                double temp_err;
                dot > 0 ? temp_err = acos(dot) : temp_err = M_PI - acos(dot);

                if(temp_err < TranslationMSEThresh){
                    inlier_new++;
                    std::vector<double> b_tilde_vec = {b_tilde(0),b_tilde(1), b_tilde(2), b_tilde(3)};
                    std::vector<double> curr_q2_fixed_vec = {q2_fixed[id].at(j)(0), q2_fixed[id].at(j)(1), q2_fixed[id].at(j)(2), q2_fixed[id].at(j)(3)};

                    ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<LineDisplacementFunctor,4,3>(
                    new LineDisplacementFunctor(b_tilde_vec, curr_q2_fixed_vec, n_vec,x_vec, R_opt));  

                    problem.AddResidualBlock(cost_function, nullptr, initial_pose.data());
                } 
            }
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
        for(int i = 0; i < plane_num; i++)
        {	
            int id = plane_id.at(i);
            int line_num = lines[id].size();
            
            Eigen::Vector4d c_new_tilde= (planes.at(id).normal_.dot(planes.at(id).x1_-t_opt) * R_opt * planes.at(id).normal_).homogeneous().normalized();

            for(int j = 0; j<line_num; j++){
                Eigen::Vector4d b_tilde = lines[id].at(j).b_.homogeneous().normalized();      
                Eigen::Vector4d q2_new = (q2_fixed[id].at(j) + b_tilde.dot(c_new_tilde) * c_new_tilde).normalized();                    
                double dot = b_tilde.dot(q2_new);
                double temp_err;
                dot > 0 ? temp_err = acos(dot) : temp_err = M_PI - acos(dot);
 
                if(temp_err < TranslationMSEThresh) current_Inlier++;
            }
        }


    if(current_Inlier > optInlier) return current_Inlier;
    else return optInlier;

    }
    
}
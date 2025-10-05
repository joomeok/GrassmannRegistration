#include "refinement.hpp"

L2LFunctor::L2LFunctor(const std::vector<double> & model_l,const std::vector<double> & data_l){
    Eigen::Vector3d model_sp = {model_l.at(0), model_l.at(1), model_l.at(2)};
    Eigen::Vector3d model_ep = {model_l.at(3), model_l.at(4), model_l.at(5)};
    Eigen::Vector3d data_sp = {data_l.at(0), data_l.at(1), data_l.at(2)};
    Eigen::Vector3d data_ep = {data_l.at(3), data_l.at(4), data_l.at(5)};
    model_d_ = model_ep - model_sp;
    Eigen::Vector3d model_m_ = model_sp.cross(model_ep);
    model_m_ = model_m_ / model_d_.norm();
    model_d_.normalize();
    model_b_ = model_d_.cross(model_m_);

    data_d_ = data_ep - data_sp;
    Eigen::Vector3d data_m_ = data_sp.cross(data_ep);
    data_m_ = data_m_ / data_d_.norm();
    data_d_.normalize();
    data_b_ = data_d_.cross(data_m_);
}

template <typename T>
bool L2LFunctor::operator()(const T* const rot, const T* const t, T* residual) const{
    Eigen::Quaternion<T> q(rot[0], rot[1], rot[2], rot[3]);
    q.normalize();
    Eigen::Matrix<T, 3, 1> trans(t[0], t[1], t[2]);
    Eigen::Matrix<T, 3, 1> data_d_jet, data_b_jet, model_d_jet, model_b_jet;
    data_d_jet << T(data_d_(0)), T(data_d_(1)), T(data_d_(2));
    data_b_jet << T(data_b_(0)), T(data_b_(1)), T(data_b_(2));
    model_d_jet << T(model_d_(0)), T(model_d_(1)), T(model_d_(2));
    model_b_jet << T(model_b_(0)), T(model_b_(1)), T(model_b_(2));
    
    Eigen::Matrix<T,3,1> d_data_R = q * data_d_jet;
    Eigen::Matrix<T,3,3> identity_jet; 
    identity_jet << T(1) , T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1);
    Eigen::Matrix<T,3,1> b_data_R_t = q * data_b_jet + q * (identity_jet - data_d_jet * data_d_jet.transpose()) * q.conjugate() * trans;

    T d_data_R_dot_d_model = d_data_R.dot(model_d_jet);
    Eigen::Matrix<T,3,1> projected_data_d = d_data_R_dot_d_model * model_d_jet;
    Eigen::Matrix<T,3,1> dir_error = projected_data_d - d_data_R;

    Eigen::Matrix<T,4,1> b_data_R_t_tilde, b_model_tilde, d_model_bar;
    b_data_R_t_tilde << b_data_R_t(0), b_data_R_t(1), b_data_R_t(2), T(1);
    b_data_R_t_tilde.normalize();  
    b_model_tilde << model_b_jet(0), model_b_jet(1), model_b_jet(2), T(1);
    b_model_tilde.normalize();
    d_model_bar << model_d_jet(0), model_d_jet(1), model_d_jet(2), T(0);
    T b_data_tilde_dot_d_model_bar = d_model_bar.dot(b_data_R_t_tilde);
    T b_data_tilde_dot_b_model_tilde = b_model_tilde.dot(b_data_R_t_tilde);
    Eigen::Matrix<T,4,1> projected_data_b_tilde = b_data_tilde_dot_d_model_bar * d_model_bar + b_data_tilde_dot_b_model_tilde * b_model_tilde;
    Eigen::Matrix<T,4,1> disp_err = projected_data_b_tilde - b_data_R_t_tilde;

    residual[0] = dir_error[0];
    residual[1] = dir_error[1];
    residual[2] = dir_error[2];
    residual[3] = disp_err[0];
    residual[4] = disp_err[1];
    residual[5] = disp_err[2];
    residual[6] = disp_err[3];
    
    return true;

}

L2PFunctor::L2PFunctor(const std::vector<double> & line,const std::vector<double> & plane){

    Eigen::Vector3d line_sp = {line.at(0), line.at(1), line.at(2)};
    Eigen::Vector3d line_ep = {line.at(3), line.at(4), line.at(5)};
    line_d_ = line_ep - line_sp;
    Eigen::Vector3d line_m_ = line_sp.cross(line_ep);
    line_m_ = line_m_ / line_d_.norm();
    line_d_.normalize();
    line_b_ = line_d_.cross(line_m_);
    
    Eigen::Vector3d x1_ = {plane.at(0), plane.at(1), plane.at(2)};
    Eigen::Vector3d x2_ = {plane.at(3), plane.at(4), plane.at(5)};
    Eigen::Vector3d x3_ = {plane.at(6), plane.at(7), plane.at(8)};

    Eigen::Vector3d normal_ = (x1_-x2_).cross(x1_-x3_).normalized();
    // Ensure normal is directed as origin -> plane
    if(normal_.dot(x1_) < 0) normal_ = -normal_;
	// Set x1-x2 as the first basis
	plane_u_ = (x1_-x2_).normalized();	
	plane_v_ = normal_.cross(plane_u_).normalized();
	plane_c_ = normal_.dot(x1_) * normal_;
}

template <typename T>
bool L2PFunctor::operator()(const T* const rot, const T* const t, T* residual) const{
    Eigen::Quaternion<T> q(rot[0], rot[1], rot[2], rot[3]);
    q.normalize();
    Eigen::Matrix<T, 3, 1> trans(t[0], t[1], t[2]);
    Eigen::Matrix<T, 3, 1> u_jet, v_jet, c_jet, d_jet, u_R, v_R, c_R_t;
    u_jet << T(plane_u_(0)), T(plane_u_(1)), T(plane_u_(2));
    v_jet << T(plane_v_(0)), T(plane_v_(1)), T(plane_v_(2));
    d_jet << T(line_d_(0)), T(line_d_(1)), T(line_d_(2));
    u_R =  q * u_jet;
    v_R = q * v_jet;
    c_jet << T(plane_c_(0)), T(plane_c_(1)), T(plane_c_(2));

    Eigen::Matrix<T,3,3> identity_jet; 
    identity_jet << T(1) , T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1);

    Eigen::Matrix<T,3,2> A; 
    A.block(0,0,3,1) = u_jet;
    A.block(0,1,3,1) = v_jet;
    c_R_t = q * c_jet + q * (identity_jet - A * A.transpose()) * q.conjugate() * trans;

    T d_dot_u_R = d_jet(0) * u_R(0) + d_jet(1) * u_R(1) + d_jet(2) * u_R(2);
    T d_dot_v_R = d_jet(0) * v_R(0) + d_jet(1) * v_R(1) + d_jet(2) * v_R(2);
    Eigen::Matrix<T,3,1> projected_d = d_dot_u_R * u_R + d_dot_v_R * v_R;
    Eigen::Matrix<T,3,1> dir_err = projected_d - d_jet;

    Eigen::Matrix<T,4,1> c_R_t_tilde, b_tilde, u_R_bar, v_R_bar;
    u_R_bar << u_R(0), u_R(1), u_R(2), T(0);
    v_R_bar << v_R(0), v_R(1), v_R(2), T(0);
    c_R_t_tilde << c_R_t(0), c_R_t(1), c_R_t(2), T(1);
    c_R_t_tilde.normalize();
    b_tilde << T(line_b_(0)), T(line_b_(1)), T(line_b_(2)), T(1);
    b_tilde.normalize();

    T u_bar_dot_b_tilde = b_tilde.dot(u_R_bar);
    T v_bar_dot_b_tilde = b_tilde.dot(v_R_bar);
    T c_R_t_tilde_dot_b_tilde = b_tilde.dot(c_R_t_tilde);
    Eigen::Matrix<T,4,1> projected_b_tilde = u_bar_dot_b_tilde * u_R_bar + v_bar_dot_b_tilde * v_R_bar + c_R_t_tilde_dot_b_tilde * c_R_t_tilde;
    Eigen::Matrix<T,4,1> disp_err = projected_b_tilde - b_tilde;

    residual[0] = dir_err[0];
    residual[1] = dir_err[1];
    residual[2] = dir_err[2];
    residual[3] = disp_err[0];
    residual[4] = disp_err[1];
    residual[5] = disp_err[2];
    residual[6] = disp_err[3];
    
    return true;
}

P2PFunctor::P2PFunctor(const std::vector<double> & model_p, const std::vector<double> & data_p){
    Eigen::Vector3d model_x1_ = {model_p.at(0), model_p.at(1), model_p.at(2)};
    Eigen::Vector3d model_x2_ = {model_p.at(3), model_p.at(4), model_p.at(5)};
    Eigen::Vector3d model_x3_ = {model_p.at(6), model_p.at(7), model_p.at(8)};

    Eigen::Vector3d model_normal_ = (model_x1_- model_x2_).cross(model_x1_- model_x3_).normalized();
    // Ensure normal is directed as origin -> plane
    if(model_normal_.dot(model_x1_) < 0) model_normal_ = -model_normal_;
	// Set x1-x2 as the first basis
	model_u_ = (model_x1_-model_x2_).normalized();	
	model_v_ = model_normal_.cross(model_u_).normalized();
	model_c_ = model_normal_.dot(model_x1_) * model_normal_;

    Eigen::Vector3d data_x1_ = {data_p.at(0), data_p.at(1), data_p.at(2)};
    Eigen::Vector3d data_x2_ = {data_p.at(3), data_p.at(4), data_p.at(5)};
    Eigen::Vector3d data_x3_ = {data_p.at(6), data_p.at(7), data_p.at(8)};

    Eigen::Vector3d data_normal_ = (data_x1_- data_x2_).cross(data_x1_- data_x3_).normalized();
    // Ensure normal is directed as origin -> plane
    if(data_normal_.dot(data_x1_) < 0) data_normal_ = -data_normal_;
	// Set x1-x2 as the first basis
	data_u_ = (data_x1_-data_x2_).normalized();	
	data_v_ = data_normal_.cross(data_u_).normalized();
	data_c_ = data_normal_.dot(data_x1_) * data_normal_;


}

template <typename T>
bool P2PFunctor::operator()(const T* const rot, const T* const t, T* residual) const{

    Eigen::Quaternion<T> q(rot[0], rot[1], rot[2], rot[3]);
    q.normalize();
    Eigen::Matrix<T, 3, 1> trans(t[0], t[1], t[2]);
    Eigen::Matrix<T, 3, 1> data_u_jet, data_v_jet, data_c_jet, model_u_jet, model_v_jet, model_c_jet, data_u_R, data_v_R, data_c_R_t;
    data_u_jet << T(data_u_(0)), T(data_u_(1)), T(data_u_(2));
    data_v_jet << T(data_v_(0)), T(data_v_(1)), T(data_v_(2));
    data_c_jet << T(data_c_(0)), T(data_c_(1)), T(data_c_(2));
    
    model_u_jet << T(model_u_(0)), T(model_u_(1)), T(model_u_(2));
    model_v_jet << T(model_v_(0)), T(model_v_(1)), T(model_v_(2));
    model_c_jet << T(model_c_(0)), T(model_c_(1)), T(model_c_(2));

    data_u_R =  q * data_u_jet;
    data_v_R = q * data_v_jet;

    Eigen::Matrix<T,3,2> A; 
    A.block(0,0,3,1) = data_u_jet;
    A.block(0,1,3,1) = data_v_jet;
    Eigen::Matrix<T,3,3> identity_jet; 
    identity_jet << T(1) , T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1);
    data_c_R_t = q * data_c_jet + q * (identity_jet - A * A.transpose()) * q.conjugate() * trans;

    T model_u_dot_data_u_R = model_u_jet.dot(data_u_R);
    T model_u_dot_data_v_R = model_u_jet.dot(data_v_R);
    Eigen::Matrix<T,3,1> projected_model_u = model_u_dot_data_u_R * data_u_R + model_u_dot_data_v_R  * data_v_R;
    Eigen::Matrix<T,3,1> dir_err1 = projected_model_u - model_u_jet;

    T model_v_dot_data_u_R = model_v_jet.dot(data_u_R);
    T model_v_dot_data_v_R = model_v_jet.dot(data_v_R);
    Eigen::Matrix<T,3,1> projected_model_v = model_v_dot_data_u_R * data_u_R + model_v_dot_data_v_R  * data_v_R;
    Eigen::Matrix<T,3,1> dir_err2 = projected_model_v - model_v_jet;

    Eigen::Matrix<T,4,1> data_c_R_t_tilde, model_c_tilde, u_R_bar, v_R_bar;
    u_R_bar << data_u_R(0), data_u_R(1), data_u_R(2), T(0);
    v_R_bar << data_v_R(0), data_v_R(1), data_v_R(2), T(0);
    data_c_R_t_tilde << data_c_R_t(0), data_c_R_t(1), data_c_R_t(2), T(1);
    data_c_R_t_tilde.normalize(); 
    model_c_tilde << model_c_jet(0), model_c_jet(1), model_c_jet(2), T(1);
    model_c_tilde.normalize();
    T model_c_dot_data_u_R = model_c_tilde.dot(u_R_bar);
    T model_c_dot_data_v_R = model_c_tilde.dot(v_R_bar);
    T model_c_dot_data_c_R_t = model_c_tilde.dot(data_c_R_t_tilde);
    Eigen::Matrix<T,4,1> projected_model_c = model_c_dot_data_u_R  * u_R_bar + model_c_dot_data_v_R * v_R_bar + model_c_dot_data_c_R_t * data_c_R_t_tilde;
    Eigen::Matrix<T,3,1> disp_err = data_c_R_t - model_c_jet;

    residual[0] = dir_err1[0];
    residual[1] = dir_err1[1];
    residual[2] = dir_err1[2];
    residual[3] = dir_err2[0];
    residual[4] = dir_err2[1];
    residual[5] = dir_err2[2];
    residual[6] = disp_err[0];
    residual[7] = disp_err[1];
    residual[8] = disp_err[2];

    return true;
}
// lines contain endpoints
std::pair<Eigen::Matrix3d, Eigen::Vector3d> L2Loptimize(std::vector<std::vector<double>> lines1, const std::vector<std::vector<double>> & lines2, const Eigen::Matrix<double,4,4>& T_init){

    int pair_num = lines1.size();
    Eigen::Matrix3d R_init = T_init.block(0,0,3,3);
    Eigen::Vector3d t_init = T_init.block(0,3,3,1);
    Eigen::Quaterniond q_init(R_init); // Construct quaternion
    Eigen::Matrix<double, 4, 1> initial_q;
    Eigen::Matrix<double, 3, 1> initial_t;
    initial_q << q_init.w(), q_init.x(), q_init.y(), q_init.z();
    initial_t << t_init(0), t_init(1), t_init(2);

    ceres::Problem problem;
    ceres::Manifold* q_manifold = new ceres::EigenQuaternionManifold();

    problem.AddParameterBlock(initial_q.data(), 4, new ceres::QuaternionManifold());
    problem.AddParameterBlock(initial_t.data(), 3);
    for (int i = 0 ; i < pair_num; i++){
        ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<L2LFunctor,7,4,3>(
        new L2LFunctor(lines1.at(i),lines2.at(i)));
        ceres::LossFunction* loss_function = new ceres::CauchyLoss(1e-3); 
        problem.AddResidualBlock(cost_function, loss_function, initial_q.data(), initial_t.data());
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    Eigen::Quaterniond q_opt(initial_q[0], initial_q[1], initial_q[2], initial_q[3]);
    Eigen::Matrix3d R_opt = q_opt.normalized().toRotationMatrix();
    Eigen::Vector3d t_opt = {initial_t[0], initial_t[1], initial_t[2]};
    // std::cout << summary.BriefReport() << std::endl;
    return std::make_pair(R_opt, t_opt);
}

std::pair<Eigen::Matrix3d, Eigen::Vector3d> L2Poptimize(std::vector<std::vector<double>> lines, const std::vector<std::vector<double>> & planes, const Eigen::Matrix<double,4,4> & T_init){

    int pair_num = lines.size();
    Eigen::Matrix3d R_init = T_init.block(0,0,3,3);
    Eigen::Vector3d t_init = T_init.block(0,3,3,1);
    Eigen::Quaterniond q_init(R_init); // Construct quaternion
    Eigen::Matrix<double, 4, 1> initial_q;
    Eigen::Matrix<double, 3, 1> initial_t;
    initial_q << q_init.w(), q_init.x(), q_init.y(), q_init.z();
    initial_t << t_init(0), t_init(1), t_init(2);

    ceres::Problem problem;
    ceres::Manifold* q_manifold = new ceres::EigenQuaternionManifold();

    problem.AddParameterBlock(initial_q.data(), 4, new ceres::QuaternionManifold());
    problem.AddParameterBlock(initial_t.data(), 3);
    for (int i = 0 ; i < pair_num; i++){
        ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<L2PFunctor,7,4,3>(
        new L2PFunctor(lines.at(i),planes.at(i)));
        ceres::LossFunction* loss_function = new ceres::CauchyLoss(1e-3); 
        problem.AddResidualBlock(cost_function, loss_function, initial_q.data(), initial_t.data());
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    Eigen::Quaterniond q_opt(initial_q[0], initial_q[1], initial_q[2], initial_q[3]);
    Eigen::Matrix3d R_opt = q_opt.normalized().toRotationMatrix();
    Eigen::Vector3d t_opt = {initial_t[0], initial_t[1], initial_t[2]};
    // std::cout << summary.BriefReport() << std::endl;
    return std::make_pair(R_opt, t_opt);
}

std::pair<Eigen::Matrix3d, Eigen::Vector3d> P2Poptimize(std::vector<std::vector<double>> planes1, const std::vector<std::vector<double>> & planes2, const Eigen::Matrix<double,4,4> & T_init){
    int pair_num = planes1.size();
    Eigen::Matrix3d R_init = T_init.block(0,0,3,3);
    Eigen::Vector3d t_init = T_init.block(0,3,3,1);
    Eigen::Quaterniond q_init(R_init); // Construct quaternion
    Eigen::Matrix<double, 4, 1> initial_q;
    Eigen::Matrix<double, 3, 1> initial_t;
    initial_q << q_init.w(), q_init.x(), q_init.y(), q_init.z();
    initial_t << t_init(0), t_init(1), t_init(2);
    // std::cout << initial_t << std::endl;
    ceres::Problem problem;
    ceres::Manifold* q_manifold = new ceres::EigenQuaternionManifold();

    problem.AddParameterBlock(initial_q.data(), 4, new ceres::QuaternionManifold());
    problem.AddParameterBlock(initial_t.data(), 3);
    for (int i = 0 ; i < pair_num; i++){
        ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<P2PFunctor,9,4,3>(
        new P2PFunctor(planes1.at(i), planes2.at(i)));
        ceres::LossFunction* loss_function = new ceres::CauchyLoss(1e-3); 
        problem.AddResidualBlock(cost_function, loss_function, initial_q.data(), initial_t.data());
    }
    
    // problem.SetParameterBlockConstant(initial_q.data());

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    // std::cout << initial_t << std::endl;

    Eigen::Quaterniond q_opt(initial_q[0], initial_q[1], initial_q[2], initial_q[3]);
    Eigen::Matrix3d R_opt = q_opt.normalized().toRotationMatrix();
    Eigen::Vector3d t_opt = {initial_t[0], initial_t[1], initial_t[2]};
    std::cout << summary.BriefReport() << std::endl;
    return std::make_pair(R_opt, t_opt);
}



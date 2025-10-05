#include "line_3d.hpp"


LINE3D::LINE3D(){
	m_ = Eigen::Vector3d::Zero();
	d_ = Eigen::Vector3d::Zero();
	b_ = Eigen::Vector3d::Zero();
}

LINE3D::LINE3D(const Eigen::Vector3d & sp, const Eigen::Vector3d & ep) : sp_(sp), ep_(ep) {
	m_ = cross(sp,ep);
	d_ = ep - sp;
	m_ = m_ / d_.norm();
	d_ = d_ / d_.norm();
	b_ = cross(d_,m_);
	// b_ = d_.cross(m_);
	empty_ = false;
}


LINE3D::LINE3D(const Eigen::Vector3d & sp, const Eigen::Vector3d & ep, const int & id) : sp_(sp), ep_(ep), plane_id(id) {
	m_ = cross(sp,ep);
	d_ = ep - sp;
	m_ = m_ / d_.norm();
	d_ = d_ / d_.norm();
	b_ = cross(d_,m_);
	// b_ = d_.cross(m_);
	empty_ = false;
}

// double ComputeLineDistance(LINE3D l1, LINE3D l2){
// 	double dot1, dot2, p_ang1, p_ang2;
// 	dot1 = l1.d_.dot(l2.d_);
// 	if(dot1 > 0) p_ang1 = std::acos(dot1);
// 	else p_ang1 = M_PI - std::acos(dot1);
// 	Eigen::Vector4d b_tilde, c_tilde;
// 	b_tilde << l1.b_(0), l1.b_(1), l1.b_(2), 1;
// 	b_tilde.normalize();
// 	c_tilde << l2.b_(0), l2.b_(1), l2.b_(2), 1;
// 	c_tilde.normalize();
// 	dot2 = b_tilde.dot(c_tilde);
// 	if(dot2 > 0) p_ang2 = std::acos(dot2);
// 	else p_ang2 = M_PI - std::acos(dot2);

// 	return pow(p_ang1,2) + pow(p_ang2,2);
// }
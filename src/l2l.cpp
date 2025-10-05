#include <stdlib.h>
#include <math.h>
#include <chrono>
#include <stdio.h>
#include <vector>
#include "l2l.hpp"
#include "line_fitting.hpp"
using namespace pmc;

GoL2L::GoL2L()
{        

	initNodeRot.a = -PI;
	initNodeRot.b = -PI;
	initNodeRot.c = -PI;
	initNodeRot.w = 2*PI;
	
	initNodeRot.l = 0;
	initNodeRot.ub = Nd;
	initNodeTrans.lb = 0;
	initNodeTrans_MC.ub = Nd;
}


Eigen::Matrix3d GoL2L::ToSkewSymmetric(const Eigen::Vector3d& rotation_vector) {
  return (Eigen::Matrix3d() << 0,                 -rotation_vector(2), rotation_vector(1),
                               rotation_vector(2), 0,                 -rotation_vector(0),
                              -rotation_vector(1), rotation_vector(0), 0).finished();
}

double GoL2L::ComputeTightRotBound(const ROTNODE & rotation_node, const Eigen::Vector3d& v, const std::vector<Eigen::Vector3d> & vertex_rotation_vectors){
	
	
	Eigen::Vector3d v_r_0 = rotation_node.SO3 * v;
    double min_cos_centre_to_vertex_angle = 1.0f;
	for (int i = 0; i < 8; i++) { // Loop over the vertices to find maximum angle
		Eigen::AngleAxis<double> angle_axis_i(vertex_rotation_vectors[i].norm(), vertex_rotation_vectors[i].normalized());
		Eigen::Matrix3d R_i = angle_axis_i.toRotationMatrix();
        double cos_centre_to_vertex_angle = v_r_0.dot(R_i *v);
      	if (cos_centre_to_vertex_angle < min_cos_centre_to_vertex_angle) min_cos_centre_to_vertex_angle = cos_centre_to_vertex_angle;
    }

	std::vector<int> edge_indices[8];
	edge_indices[0] = {1, 2, 4}; // Vertex 0 connects to vertices 1, 2 and 4
	edge_indices[1] = {3, 5};
	edge_indices[2] = {3, 6};
	edge_indices[3] = {7};
	edge_indices[4] = {5, 6};
	edge_indices[5] = {7};
	edge_indices[6] = {7};
	edge_indices[7] = {};
	double rotation_uncertainty_radius = maxRotAng[rotation_node.l];

	for(int i = 0; i< 8; i++){
		Eigen::Vector3d r_i = vertex_rotation_vectors[i];
		Eigen::AngleAxis<double> angle_axis_i(r_i.norm(), r_i.normalized());
		Eigen::Matrix3d R_i = angle_axis_i.toRotationMatrix();
		Eigen::Vector3d v_r_i = R_i * v;

		for (std::vector<int>::iterator it = edge_indices[i].begin(); it != edge_indices[i].end(); ++it) {
			Eigen::Vector3d r_j = vertex_rotation_vectors[*it];
			Eigen::AngleAxis<double> angle_axis_j(r_j.norm(), r_j.normalized());
			Eigen::Matrix3d R_j = angle_axis_j.toRotationMatrix();
			Eigen::Vector3d v_r_j = R_j * v;

			double derivative_at_0 = v.transpose() * rotation_node.SO3.transpose() * R_i * ToSkewSymmetric(v) * (r_i * r_i.transpose() + (R_i.transpose() - Eigen::Matrix3d::Identity()) * ToSkewSymmetric(r_i)) * (r_j - r_i);
			if(derivative_at_0 >= 0.0f){
				double derivative_at_1 = v.transpose() * rotation_node.SO3.transpose() * R_j * ToSkewSymmetric(v) * (r_j * r_j.transpose() + (R_j.transpose() - Eigen::Matrix3d::Identity()) * ToSkewSymmetric(r_j)) * (r_j - r_i);
				if(derivative_at_1 <=0.0f){
					
				static const float inverted_golden_ratio = 0.5f * (std::sqrt(5) - 1.0f);
				static const float one_minus_inverted_golden_ratio = 0.5f * (3.0f - std::sqrt(5));
				static const float tol = M_PI / 2048.0f; // at most 0.088 degrees incorrect
				Eigen::Vector3d a = r_i;
				Eigen::Vector3d b = r_j;
				double fa = v_r_0.dot(v_r_i);
				double fb = v_r_0.dot(v_r_j); 
				Eigen::Vector3d c = b - inverted_golden_ratio * (b - a);
				Eigen::Vector3d d = a + inverted_golden_ratio * (b - a);
				
				Eigen::AngleAxis<double> angle_axis_c(c.norm(), c.normalized());
				Eigen::Matrix3d R_c = angle_axis_c.toRotationMatrix();
				
				Eigen::AngleAxis<double> angle_axis_d(d.norm(), d.normalized());
				Eigen::Matrix3d R_d = angle_axis_d.toRotationMatrix();

				double fc = v_r_0.dot(R_c * v);
				double fd = v_r_0.dot(R_d * v);

				while ((c - d).norm() > tol) {
					if (fc < fd) {
						b = d;
						fb = fd;
						d = c;
						fd = fc;
						c = a + one_minus_inverted_golden_ratio * (b - a);
						Eigen::AngleAxis<double> angle_axis_c(c.norm(), c.normalized());
						Eigen::Matrix3d R_c = angle_axis_c.toRotationMatrix();
						fc = v_r_0.dot(R_c * v);
					} else {
						a = c;
						fa = fc;
						c = d;
						fc = fd;
						d = b - one_minus_inverted_golden_ratio * (b - a);

						Eigen::AngleAxis<double> angle_axis_d(d.norm(), d.normalized());
						Eigen::Matrix3d R_d = angle_axis_d.toRotationMatrix();
						fd = v_r_0.dot(R_d * v);
					}
				}
				c = 0.5f * (b + a);
				Eigen::AngleAxis<double> angle_axis_temp(c.norm(), c.normalized());
				Eigen::Matrix3d R_temp = angle_axis_temp.toRotationMatrix();
				fc = v_r_0.dot(R_temp * v);
				if (fc < min_cos_centre_to_vertex_angle) min_cos_centre_to_vertex_angle = fc;
				}
			}
		}
	}

	return std::min(std::acos(min_cos_centre_to_vertex_angle), rotation_uncertainty_radius);
}





void GoL2L::Initialize()
{
	int i, j;
	double maxAngle;

	for(i = 0; i < 30; i++)
	{
		double sigma = initNodeRot.w/pow(2.0,i)/2.0; // Half-side length of each level of rotation subcube
		maxAngle = SQRT3*sigma;
		if(maxAngle > PI)
			maxAngle = PI;
		maxRotAng.push_back(maxAngle);
	}

	// Temporary Variables
	psi1 = new double[Nd];
	is_positive_lb = new bool[Nd];

	initNodeRot.ub = Nd;
	optNodeRot = initNodeRot;
	optNodeTrans = initNodeTrans;
	optNodeTrans_MC = initNodeTrans_MC;

	// Initialise so-far-best rotation and translation matrices
	optR = Eigen::Matrix3d::Identity();
	optT = Eigen::Vector3d::Zero();

	inlierNum = Nd;

	RotationSSEThresh = RotationMSEThresh * inlierNum;
}

void GoL2L::Clear()
{
	delete(psi1);
	delete(is_positive_lb);
}

void GoL2L::TranslationSearch(std::vector<std::pair<int,int>> Correspondence, std::vector<int> MaxClique){
	std::vector<LINE3D> inlier_data, inlier_model;
	int c_num = MaxClique.size();
	
	std::cout << 0 << std::endl;
	
	for(const auto & idx : MaxClique){
		inlier_data.push_back(lData.at(Correspondence.at(idx).first)); 
		inlier_model.push_back(lModel.at(Correspondence.at(idx).second));
	}

	std::cout << 1 << std::endl;

	float transX, transY, transZ;
    float lb, ub, optErrorT, error;
	TRANSNODE nodeTrans, nodeTransParent;
    std::priority_queue<TRANSNODE> queueTrans;
	queueTrans.push(initNodeTrans);
	long long count = 0;
	Eigen::Vector3d t_init;
	optErrorT = optimizer_.optimize_t(inlier_data, inlier_model, optT, t_init, optR);
	optT = t_init;
	TranslationSSEThresh = TranslationMSEThresh * inlier_data.size();
	while(1)
	{
		if(queueTrans.empty()){
			break;
		}

		nodeTransParent = queueTrans.top();
		queueTrans.pop();
		if(optErrorT-nodeTransParent.lb < TranslationSSEThresh)
		{	
			break;
		}
			count++;

        nodeTrans.w = nodeTransParent.w/2;
        for(int j = 0; j < 8; j++)
		{   
            std::vector<double> current_distance;
            nodeTrans.x = nodeTransParent.x + (j&1)*nodeTrans.w ;
			nodeTrans.y = nodeTransParent.y + (j>>1&1)*nodeTrans.w ;
			nodeTrans.z = nodeTransParent.z + (j>>2&1)*nodeTrans.w ;
            
            transX = nodeTrans.x + nodeTrans.w/2;
			transY = nodeTrans.y + nodeTrans.w/2;
			transZ = nodeTrans.z + nodeTrans.w/2;

            Eigen::Vector3d t0 = {transX, transY, transZ};
            std::vector<double> maxTransDist;
            for(int i = 0; i< c_num; i++)
			{
                Eigen::Vector3d b_i_t0_vec = inlier_data.at(i).b_ - inlier_data.at(i).d_.cross(t0.cross(inlier_data.at(i).d_));
                Eigen::Vector3d b_i_R_t0_vec = optR * b_i_t0_vec;
                current_distance.push_back(ComputeDispDist(b_i_R_t0_vec, inlier_model.at(i).b_));
				
                double max_trans_dist = 0;
                for(int k =0; k < 8; k++){
                    Eigen::Vector3d vertice;
                    vertice << nodeTrans.x + (k&1)*nodeTrans.w, nodeTrans.y + (k>>1&1)*nodeTrans.w, nodeTrans.z + (k>>2&1)*nodeTrans.w;
                    Eigen::Vector3d b_i_R_vert = optR * (lData.at(i).b_ - lData.at(i).d_.cross(vertice.cross(lData.at(i).d_)));
					double dist = std::acos((b_i_R_t0_vec.dot(b_i_R_vert) + 1) / (sqrt(1 + pow(b_i_R_t0_vec.norm(),2)) * sqrt(1 + pow(b_i_R_vert.norm(),2))));
					if(dist > max_trans_dist) max_trans_dist = dist;
                }
                if(max_trans_dist > M_PI /2) max_trans_dist = M_PI / 2;
				maxTransDist.push_back(max_trans_dist);	     
			}

            ub = 0;
            for(int i = 0; i < c_num; i++)
            {	
                ub += current_distance.at(i) * current_distance.at(i);
            }  
			if(ub < optErrorT)
			{
				// Update optimal error and rotation/translation nodes
				optNodeTrans = nodeTrans;


				// Run ICP
				t_init << optNodeTrans.x+optNodeTrans.w/2, optNodeTrans.y+optNodeTrans.w/2, optNodeTrans.z+optNodeTrans.w/2;
				Eigen::Vector3d t_optimized;
				error = optimizer_.optimize_t(inlier_data,inlier_model,t_init,t_optimized,optR);
				
				if(error < optErrorT)
				{	
					optErrorT = error;
					optT = t_optimized;
				}


				std::priority_queue<TRANSNODE> queueTransNew;

		        while(!queueTrans.empty())
				{
					TRANSNODE node = queueTrans.top();
					queueTrans.pop();
					if(node.lb < optErrorT)
						queueTransNew.push(node);
					else
						break;
				}
				queueTrans = queueTransNew;
			}

            lb = 0;
            for(int i = 0; i < c_num; i++){
                double psi_temp = current_distance.at(i) - maxTransDist[i];
                if(psi_temp < 0) psi_temp = 0;
                lb += psi_temp * psi_temp;
            }

            if(lb >= optErrorT)
			{
				//discard
				continue;
			}

			nodeTrans.ub = ub;
			nodeTrans.lb = lb;
			queueTrans.push(nodeTrans);
		}
		if(count > 1000) break;

	}
}



std::vector<std::pair<int,int>> GoL2L::RotationSearch(){
	ROTNODE nodeRot, nodeRotParent;
	float lb, ub, error;
	clock_t clockBeginICL;
	

	// std::priority_queue<RotNode> queueRot;
	std::priority_queue<ROTNODE> queueRot;
	float v1, v2, v3, t, ct, ct2,st, st2;
	float tmp121, tmp122, tmp131, tmp132, tmp231, tmp232;
	float R11, R12, R13, R21, R22, R23, R31, R32, R33;

	optInlier = -1;
	std::vector<std::pair<int,int>> correspondence;

	Eigen::Matrix3d R_init;
	optimizer_.optimize_R(lData, lModel, optR, R_init, RotationMSEThresh, optInlier);
	
	int current_Inlier = 0;
	for(int i = 0; i < Nd; i++)
	{	
		std::vector<double> data_direction;
		Eigen::Vector3d lData_R = optR * lData[i].d_;
		Eigen::Vector3d lModel_R = lModel[i].d_;
		double dot = lData_R.dot(lModel_R);
		double error;
		dot > 0 ? error = acos(dot) : error = PI - acos(dot);
		if(error < RotationMSEThresh) current_Inlier++;
	}
	optInlier = current_Inlier;
	


	queueRot.push(initNodeRot);
	long long count = 0;

	while(1)
	{	
		if(queueRot.empty())
		{
		  std::cout << "Rotation Queue Empty" << std::endl;
		  std::cout << "Inlier num: " << optInlier << ", UB: " << ub << std::endl;
		  break;
		}

		// Access rotation cube with lowest lower bound...
		nodeRotParent = queueRot.top();
		// ...and remove it from the queue
		queueRot.pop();


		if(nodeRotParent.ub <= optInlier)
		{	
			break;
		}


		count++;
		nodeRot.w = nodeRotParent.w/2;
		nodeRot.l = nodeRotParent.l+1;

		for(int j = 0; j < 8; j++)
		{	
			ub = 0;
			lb = 0;
			nodeRot.a = nodeRotParent.a + (j&1)*nodeRot.w ;
			nodeRot.b = nodeRotParent.b + (j>>1&1)*nodeRot.w ;
			nodeRot.c = nodeRotParent.c + (j>>2&1)*nodeRot.w ;
			v1 = nodeRot.a + nodeRot.w/2;
			v2 = nodeRot.b + nodeRot.w/2;
			v3 = nodeRot.c + nodeRot.w/2;
			
			if(sqrt(v1*v1+v2*v2+v3*v3)-SQRT3*nodeRot.w/2 > PI)
			{
				continue;
			}

			t = sqrt(v1*v1 + v2*v2 + v3*v3);
			if(t > 0)
			{
				v1 /= t;
				v2 /= t;
				v3 /= t;

				ct = cos(t);
				ct2 = 1 - ct;
				st = sin(t);
				st2 = 1 - st;

				tmp121 = v1*v2*ct2; tmp122 = v3*st;
				tmp131 = v1*v3*ct2; tmp132 = v2*st;
				tmp231 = v2*v3*ct2; tmp232 = v1*st;

				R11 = ct + v1*v1*ct2;		R12 = tmp121 - tmp122;		R13 = tmp131 + tmp132;
				R21 = tmp121 + tmp122;		R22 = ct + v2*v2*ct2;		R23 = tmp231 - tmp232;
				R31 = tmp131 - tmp132;		R32 = tmp231 + tmp232;		R33 = ct + v3*v3*ct2;
			}
			nodeRot.SO3 << R11, R12, R13, R21, R22, R23, R31, R32, R33;

			for(int i = 0; i<Nd ; i++){                
				Eigen::Vector3d d_i_R0_vec = nodeRot.SO3 * lData.at(i).d_;
				Eigen::Vector3d lModel_d = lModel.at(i).d_;
				double dot = d_i_R0_vec.dot(lModel_d);
				double error;
				dot > 0 ? error = acos(dot) : error = PI - acos(dot);

				psi1[i] = error;

				if(RotationMSEThresh - psi1[i] > 0){
					lb++;
				}
			}

			if(optInlier < 2*lb){
				optNodeRot = nodeRot;
				R_init = nodeRot.SO3;
				Eigen::Matrix3d R_optimized;
				int current_Inlier = optimizer_.optimize_R(lData, lModel, R_init, R_optimized, RotationMSEThresh, optInlier);


				if(optInlier < current_Inlier)
				{	
					optInlier = current_Inlier;
					optR = R_optimized;
					// Update Correspondence Pair for MC 
					std::vector<std::pair<int,int>> current_correspondence;
					for (int i=0; i< Nd; i++){
						std::vector<double> data_direction;
						Eigen::Vector3d lData_R = optR * lData[i].d_;
						Eigen::Vector3d lModel_d = lModel[i].d_;
						double dot = lData_R.dot(lModel_d);
						double error;
						dot > 0 ? error = acos(dot) : error = PI - acos(dot);

						if(error < 0.03){

							current_correspondence.push_back(std::make_pair(i,i));
						}
					}
					correspondence = current_correspondence;
				}

				std::priority_queue<ROTNODE> queueRotNew;

		        while(!queueRot.empty())
				{
					// RotNode node = queueRot.top();
					ROTNODE node = queueRot.top();
					queueRot.pop();
					if(node.ub > optInlier)
						queueRotNew.push(node);
					else
						break;
				}
				queueRot = queueRotNew;
			}



			for(int i = 0; i<Nd ; i++){    
				double current_maxRot;


				if(maxRotAng[nodeRot.l] >= M_PI/2) current_maxRot = M_PI /2;
				else current_maxRot = maxRotAng[nodeRot.l];  

				if(RotationMSEThresh - psi1[i] + current_maxRot > 0){
					ub++;
				} 
			}

			if(ub < optInlier)
			{
				continue;
			}

    		nodeRot.ub = ub;
			nodeRot.lb = lb;
			queueRot.push(nodeRot);
		}

	}
	return correspondence;
}

void GoL2L::TranslationMCSearch(std::vector<std::pair<int,int>> correspondence){
	float transX, transY, transZ;
    float lb, ub, optErrorT, error;
	TRANSNODE_MC nodeTrans, nodeTransParent;
    std::priority_queue<TRANSNODE_MC> queueTrans;
	initNodeTrans_MC.ub = Nd;
	queueTrans.push(initNodeTrans_MC);
	
	Eigen::Vector3d t_init; 
	int matched_line_num = correspondence.size();

	std::vector<std::vector<Eigen::Vector3d>> tot_candidate_model_disps;
	std::vector<std::vector<LINE3D>> tot_candidate_models;
	std::vector<LINE3D> lData_MC;
	std::vector<LINE3D> lModel_MC;

	double * psi2 = new double[matched_line_num];

	for(int i = 0; i < matched_line_num; i++){
		int inlier_idx = correspondence.at(i).first;
		lData_MC.push_back(lData.at(inlier_idx));
		lModel_MC.push_back(lModel.at(inlier_idx));

	}




	optInlier = -1;
	optimizer_.optimize_t_MC(lData_MC, lModel_MC, optT, t_init, optR, TranslationMSEThresh, optInlier);
	

    int current_Inlier = 0;
    for(int i = 0; i < matched_line_num; i++)
    {	
        Eigen::Vector3d data_disp_Rt = optR * (lData_MC.at(i).b_ - lData_MC.at(i).d_.cross(optT.cross(lData_MC.at(i).d_)));
		
		double dot = data_disp_Rt.homogeneous().normalized().dot(lModel_MC.at(i).b_.homogeneous().normalized());
		double error;
		dot > 0 ? error = acos(dot) : error = PI - acos(dot);

        if(error < TranslationMSEThresh) current_Inlier++;
    }

	optInlier = current_Inlier;

	// optInlier = -1;

	long long count = 0;

	while(1)
	{	
		if(queueTrans.empty())
		{
		  std::cout << "Translation Queue Empty" << std::endl;
		  std::cout << "Inlier num: " << optInlier << ", UB: " << ub << std::endl;
		  break;
		}

		// Access rotation cube with lowest lower bound...
		nodeTransParent = queueTrans.top();
		// ...and remove it from the queue
		queueTrans.pop();

		// Exit if the optError is less than or equal to the lower bound plus a small epsilon
		if(count > 0 && nodeTransParent.ub <= optInlier)
		{	
			break;
		}

		

		nodeTrans.w = nodeTransParent.w/2;
		// For inital division to 8 points
        for(int j = 0; j < 8; j++)
		{   

			ub = 0;
			lb = 0;
			count++;
            std::vector<double> current_distance;
            nodeTrans.x = nodeTransParent.x + (j&1)*nodeTrans.w ;
			nodeTrans.y = nodeTransParent.y + (j>>1&1)*nodeTrans.w ;
			nodeTrans.z = nodeTransParent.z + (j>>2&1)*nodeTrans.w ;
            
            transX = nodeTrans.x + nodeTrans.w/2;
			transY = nodeTrans.y + nodeTrans.w/2;
			transZ = nodeTrans.z + nodeTrans.w/2;

            Eigen::Vector3d t0 = {transX, transY, transZ};
			std::vector<double> maxTransDist;
			
			for(int i = 0; i< matched_line_num; i++){
                Eigen::Vector3d b_i_t0_vec = lData_MC.at(i).b_ - lData_MC.at(i).d_.cross(t0.cross(lData_MC.at(i).d_));
                Eigen::Vector3d b_i_R_t0_vec = optR * b_i_t0_vec;
				std::vector<double> b_i_R_t0 = {b_i_R_t0_vec(0), b_i_R_t0_vec(1),b_i_R_t0_vec(2)};
				double dot = b_i_R_t0_vec.homogeneous().normalized().dot(lModel_MC.at(i).b_.homogeneous().normalized());
				double dist;
				dot > 0 ? dist = acos(dot) : dist = PI - acos(dot);


				psi2[i] = dist;

				if(TranslationMSEThresh - psi2[i] > 0){
					lb++;
				}

                double max_trans_dist = 0;
                for(int k =0; k < 8; k++){
                    Eigen::Vector3d vertice;
                    vertice << nodeTrans.x + (k&1)*nodeTrans.w, nodeTrans.y + (k>>1&1)*nodeTrans.w, nodeTrans.z + (k>>2&1)*nodeTrans.w;
                    Eigen::Vector3d b_i_R_vert = optR * (lData_MC.at(i).b_ - lData_MC.at(i).d_.cross(vertice.cross(lData_MC.at(i).d_)));
					double dist = std::acos((b_i_R_t0_vec.dot(b_i_R_vert) + 1) / (sqrt(1 + pow(b_i_R_t0_vec.norm(),2)) * sqrt(1 + pow(b_i_R_vert.norm(),2))));
					if(dist > max_trans_dist) max_trans_dist = dist;
                }
                if(max_trans_dist > M_PI /2) max_trans_dist = M_PI / 2;
				maxTransDist.push_back(max_trans_dist);	   
			}

			if(optInlier < 2*lb){
				optNodeTrans_MC = nodeTrans;
				t_init << optNodeTrans_MC.x+optNodeTrans_MC.w/2, optNodeTrans_MC.y+optNodeTrans_MC.w/2, optNodeTrans_MC.z+optNodeTrans_MC.w/2;
				Eigen::Vector3d t_optimized;

				int current_Inlier = optimizer_.optimize_t_MC(lData_MC, lModel_MC, t_init, t_optimized, optR, TranslationMSEThresh, optInlier);

				// std::cout << "optInlier : " << optInlier << " current Inlier : " << current_Inlier << std::endl;
				if(optInlier < current_Inlier)
				{	
					optInlier = current_Inlier;
					optT = t_optimized;
					// Update Correspondence Pair for MC 
				}

				std::priority_queue<TRANSNODE_MC> queueTransNew;

		        while(!queueTrans.empty())
					{
					TRANSNODE_MC node = queueTrans.top();
					queueTrans.pop();
					if(node.ub > optInlier)
						queueTransNew.push(node);
					else
						break;
				}
				queueTrans = queueTransNew;
			}

			for(int i = 0; i< matched_line_num ; i++){     
				// std::cout << "psi2 : " <<  psi2[i] << " maxTransDist : " << maxTransDist.at(i) << std::endl;
				if(TranslationMSEThresh - psi2[i] + maxTransDist.at(i) > 0){
					ub++;
				} 
			}
			// std::cout << "UB : " << ub << "LB : " << lb << std::endl;

			if(ub < optInlier)
			{
				continue;
			}

    		nodeTrans.ub = ub;
			nodeTrans.lb = lb;
			queueTrans.push(nodeTrans);
		}


	}
	
	delete psi2;
}
void convertToMMCFFile(bool** graph, const std::string& filename, int nonZeroCount, int N) {

    // Open file
    std::ofstream file(filename);

    // Check if file is open
    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Write header
    file << "%%MatrixMarket matrix coordinate pattern symmetric" << std::endl;

    // Write dimensions and nonzero count
    file << N << " " << N << " " << nonZeroCount << std::endl;

    // Write nonzero entries
    for (int i = 0; i < N; ++i) {
        for (int j = i; j < N; ++j) {
            if (graph[i][j]) {
                file << j + 1 << " " << i + 1 << std::endl;
            }
        }
    }
    // Close file
    file.close();
}

std::vector<int> GoL2L::FindMaxClique(const std::vector<std::pair<int,int>> & correspondence){
	int pair_num = correspondence.size();
	// bool** consistency_graph = new bool [pair_num][pair_num];
	std::cout << "Pair num : " << pair_num << std::endl;
	bool** consistency_graph = new bool*[pair_num];
	for(int i = 0; i < pair_num; ++i){
		consistency_graph[i] = new bool[pair_num];
	}
	for(int i = 0; i< pair_num; i++) consistency_graph[i][i] = false;
	int nonZeroCount = 0;
	for (auto it1 = correspondence.begin(); it1 != correspondence.end(); ++it1) {
		int idx1 = std::distance(correspondence.begin(), it1);
		auto it2 = it1;
		++it2;

		Eigen::Vector3d C1_data_d = lData.at(it1->first).d_;
		Eigen::Vector3d C1_model_d = lModel.at(it1->second).d_;

		Eigen::Vector3d C1_data_P = lData.at(it1->first).sp_;
		Eigen::Vector3d C1_model_P = lModel.at(it1->second).sp_;

		for (; it2 != correspondence.end(); ++it2) {
			int idx2 = std::distance(correspondence.begin(), it2);
			Eigen::Vector3d C2_data_d = lData.at(it2->first).d_;
			Eigen::Vector3d C2_model_d = lModel.at(it2->second).d_;

			Eigen::Vector3d C2_data_P = lData.at(it2->first).sp_;
			Eigen::Vector3d C2_model_P = lModel.at(it2->second).sp_;
			double dist1 = std::abs((C2_data_P - C1_data_P).dot(C1_data_d.cross(C2_data_d))) / C1_data_d.cross(C2_data_d).norm();

			double dist2 = std::abs((C2_model_P - C1_model_P).dot(C1_model_d.cross(C2_model_d))) / C1_model_d.cross(C2_model_d).norm();
			// std::cout << std::abs(dist1 - dist2)  << std::endl;
			if(std::abs(dist1 - dist2) < 0.3){
				consistency_graph[idx1][idx2] = true; consistency_graph[idx2][idx1] = true; nonZeroCount++;
			} 
			else{
				consistency_graph[idx1][idx2] = false; 	consistency_graph[idx2][idx1] = false;
			} 


		}
	}


	convertToMMCFFile(consistency_graph, "./Graph.mtx", nonZeroCount, pair_num);

	int argc = 5;
	char *argv[argc];
	argv[1] = (char*)" ";
	argv[1] = (char*)"-f";
	argv[2] = (char*)"./Graph.mtx";
	argv[3] = (char*)"-a";
	argv[4] = (char*)"0";

	input in(argc, argv);

    //! read graph
    pmc_graph G(in.graph_stats,in.graph);
    if (in.graph_stats) { G.bound_stats(in.algorithm, in.lb, G); }

    //! ensure wait time is greater than the time to recompute the graph data structures
    if (G.num_edges() > 1000000000 && in.remove_time < 120)  in.remove_time = 120;
    else if (G.num_edges() > 250000000 && in.remove_time < 10) in.remove_time = 10;

    //! upper-bound of max clique
    double seconds = get_time();
    G.compute_cores();
    if (in.ub == 0) {
        in.ub = G.get_max_core() + 1;
    }

    //! lower-bound of max clique
    std::vector<int> C;
    if (in.lb == 0 && in.heu_strat != "0") { // skip if given as input
        pmc_heu maxclique(G,in);
        in.lb = maxclique.search(G, C);

    }

	// C contains max clique


    else if (in.algorithm >= 0) {
        switch(in.algorithm) {
            case 0: {
                //! k-core pruning, neigh-core pruning/ordering, dynamic coloring bounds/sort
                if (G.num_vertices() < in.adj_limit) {
                    G.create_adj();
                    pmcx_maxclique finder(G,in);
                    finder.search_dense(G,C);
                    break;
                }
                else {
                    pmcx_maxclique finder(G,in);
                    finder.search(G,C);
                    break;
                }
            }
            case 1: {
                //! k-core pruning, dynamic coloring bounds/sort
                if (G.num_vertices() < in.adj_limit) {
                    G.create_adj();
                    pmcx_maxclique_basic finder(G,in);
                    finder.search_dense(G,C);
                    break;
                }
                else {
                    pmcx_maxclique_basic finder(G,in);
                    finder.search(G,C);
                    break;
                }
            }
            case 2: {
                //! simple k-core pruning (four new pruning steps)
                pmc_maxclique finder(G,in);
                finder.search(G,C);
                break;
            }
            default:
                std::cout << "algorithm " << in.algorithm << " not found." << std::endl;
                break;
        }


        if (C.size() < in.param_ub)
            std::cout << "Clique of size " << in.param_ub << " does not exist." << std::endl;
    }
	
	for(int i = 0; i < pair_num; ++i) {
		delete [] consistency_graph[i];
	}
	delete [] consistency_graph;

	return C;
}

double GoL2L::ComputeDispDist(Eigen::Vector3d disp1, Eigen::Vector3d disp2){
		Eigen::Vector4d b1, b2;
		b1 << disp1(0), disp1(1), disp1(2), 1;
		b2 << disp2(0), disp2(1), disp2(2), 1;
		b1 = b1 / b1.norm();
		b2 = b2 / b2.norm();
		// double acos = std::acos(b1.dot(b2));
		double acos;
		b1.dot(b2) > 0 ? acos = std::acos(b1.dot(b2)) : acos = M_PI - std::acos(b1.dot(b2));
		double distc = acos;
		return distc;
}

double GoL2L::ComputeLineDist(LINE3D l1, LINE3D l2){
		Eigen::Vector3d d1, d2;
		d1 << l1.d_(0), l1.d_(1), l1.d_(2);
		d2 << l2.d_(0), l2.d_(1), l2.d_(2);
		double acos1;

		d1.dot(d2) > 0 ? acos1 = std::acos(d1.dot(d2)) : acos1 = M_PI - std::acos(d1.dot(d2));

		Eigen::Vector4d b1, b2;
		b1 << l1.b_(0), l1.b_(1), l1.b_(2), 1;
		b2 << l2.b_(0), l2.b_(1), l2.b_(2), 1;
		b1 = b1 / b1.norm();
		b2 = b2 / b2.norm();
		// double acos = std::acos(b1.dot(b2));
		double acos2;
		b1.dot(b2) > 0 ? acos2 = std::acos(b1.dot(b2)) : acos2 = M_PI - std::acos(b1.dot(b2));
        double distc = pow(acos1,2) + pow(acos2,2);
        return distc;		
}



std::pair<double,int> GoL2L::Register()
{	
	Initialize();

    std::chrono::high_resolution_clock::time_point t_time_begin, t_time_end;
    std::chrono::microseconds t_diff_usec;
    t_time_begin = std::chrono::high_resolution_clock::now();

	auto correspondence = RotationSearch();
	std::cout <<optR << std::endl;
	auto MaxClique = FindMaxClique(correspondence);
	if(MaxClique.size() > 2) TranslationSearch(correspondence,MaxClique);
	else TranslationMCSearch(correspondence);
	// TranslationSearch(correspondence,MaxClique);	
	// TranslationMCSearch(correspondence);

    t_time_end = std::chrono::high_resolution_clock::now();
    t_diff_usec = std::chrono::duration_cast<std::chrono::microseconds>(t_time_end - t_time_begin);             
    double t_diff_sec = t_diff_usec.count() / 1000000.0;

	printf("time : %.6fs\n",t_diff_sec);


	Clear();

	return std::make_pair(t_diff_sec * 1000, optInlier);
}

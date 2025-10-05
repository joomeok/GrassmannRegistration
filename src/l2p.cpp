#include "l2p.hpp"


#include <fstream>

GoL2P::GoL2P()
{        
	initNodeRot.a = -PI;
	initNodeRot.b = -PI;
	initNodeRot.c = -PI;
	initNodeRot.w = 2*PI;
	initNodeRot.l = 0;
	initNodeRot.ub = Nm;

	initNodeRot_full.a = -PI;
	initNodeRot_full.b = -PI;
	initNodeRot_full.c = -PI;
	initNodeRot_full.w = 2*PI;
	initNodeRot_full.lb = 0;

	initNodeTrans.lb = 0;
	
}



void GoL2P::Initialize()
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
	psi1 = new std::vector<double>[Nd];

	// Initialise so-far-best rotation and translation nodes
	initNodeRot.ub = Nm;
	optNodeRot = initNodeRot;
	optNodeRot_full = initNodeRot_full;
	
	optNodeTrans = initNodeTrans;
	optNodeTrans_MC = initNodeTrans_MC;

	// Initialise so-far-best rotation and translation matrices
	optR = Eigen::Matrix3d::Identity();
	optT = Eigen::Vector3d::Zero();

	inlierNum = Nm;
	RotationSSEThresh = RotationMSEThresh * inlierNum;
}

void GoL2P::Clear()
{
	delete [] psi1;
}


double GoL2P::RotationSearch_full(){
	ROTNODE_FULL nodeRot, nodeRotParent;
	clock_t clockBeginL2P;
	
	std::priority_queue<ROTNODE_FULL> queueRot;
	float v1, v2, v3, t, ct, ct2,st, st2;
	float tmp121, tmp122, tmp131, tmp132, tmp231, tmp232;
	float R11, R12, R13, R21, R22, R23, R31, R32, R33;
    float lb, ub, optErrorR, error;
	optErrorR = 1e3;
	queueRot.push(initNodeRot_full);


	Eigen::Matrix3d R_init;
	optErrorR = optimizer_.optimize_R_full(RotationMSEThresh, lModel, pData, optR, R_init, optErrorR, Nm);
	

	optR = R_init;
	RotationSSEThresh = RotationMSEThresh * Nm;

	long long count = 0;

	while(1)
	{	
		if(queueRot.empty())
		{
		  std::cout << "Rotation Queue Empty" << std::endl;
		  std::cout << "Inlier num: " << optInlier_ << ", UB: " << ub << std::endl;
		  break;
		}

		// Access rotation cube with lowest lower bound...
		nodeRotParent = queueRot.top();
		// ...and remove it from the queue
		queueRot.pop();
		// Exit if the optError is less than or equal to the lower bound plus a small epsilon
		if(optErrorR-nodeRotParent.lb < RotationSSEThresh)
		{	
			std::cout << "OptErrorR: " << optErrorR << " nodeRotParent.lb " << nodeRotParent.lb << " SSE " << RotationSSEThresh << std::endl;
			break;
		}

		// if(count>0 && count%300 == 0)
			// printf("LB=%f  L=%d\n",nodeRotParent.lb, nodeRotParent.l);
		
		count ++;
		nodeRot.w = nodeRotParent.w/2;
		nodeRot.l = nodeRotParent.l+1;

		for(int j = 0; j < 8; j++)
		{	

			nodeRot.a = nodeRotParent.a + (j&1)*nodeRot.w;
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
	

			ub = 0;
			for(int i = 0; i<Nd; i++){                
				int line_num = lModel[i].size();
				Eigen::Vector3d u = pData.at(i).u_;
				Eigen::Vector3d v = pData.at(i).v_;
				for(int j = 0; j<line_num; j++){
					Eigen::Vector3d d = lModel[i].at(j).d_;
					Eigen::Vector3d q_i_R0_vec = (d.dot(nodeRot.SO3 *u) * nodeRot.SO3 * u + d.dot(nodeRot.SO3 *v) * nodeRot.SO3 * v).normalized(); 
					double dot = d.dot(q_i_R0_vec);
					double dist;
					dot > 0 ? dist = safe_acos(dot) : dist = M_PI - safe_acos(dot);
					psi1[i].push_back(dist);
					ub += dist * dist; 
				}
			}

			if(ub < optErrorR){
				optNodeRot_full = nodeRot;
				R_init = nodeRot.SO3;
				Eigen::Matrix3d R_optimized;
				error = optimizer_.optimize_R_full(RotationMSEThresh, lModel, pData, R_init, R_optimized, optErrorR, Nm);

				if(error < optErrorR)
				{	
					optErrorR = error;
					optR = R_optimized;

				}

				std::priority_queue<ROTNODE_FULL> queueRotNew;

		        while(!queueRot.empty())
				{
					// RotNode node = queueRot.top();
					ROTNODE_FULL node = queueRot.top();
					queueRot.pop();
					if(node.lb < optInlier_)
						queueRotNew.push(node);
					else
						break;
				}
				queueRot = queueRotNew;
			}
			lb = 0 ;
			for(int i = 0; i<Nd ; i++){       
				double current_maxRot;
				if(maxRotAng[nodeRot.l] >= M_PI/2) current_maxRot = M_PI /2;
				else current_maxRot = maxRotAng[nodeRot.l];          
				int line_num = lModel[i].size();
				for(int j = 0; j <line_num; j++){
					double d_R0_dist = psi1[i].at(j);
					double psi_temp = psi1[i].at(j) - current_maxRot;
					if(psi_temp <0) psi_temp =0;
					lb += psi_temp * psi_temp;
				}
			}

			for(int i = 0; i<Nd ; i++) psi1[i].clear();
			if(lb >= optErrorR)
			{
				continue;
			}

    		nodeRot.ub = ub;
			nodeRot.lb = lb;
			queueRot.push(nodeRot);
		}
	}
	return optErrorR;
}




std::map<int,std::vector<LINE3D>> GoL2P::RotationSearch(){
	ROTNODE nodeRot, nodeRotParent;
	float lb, ub, error;
	clock_t clockBeginL2P;
	
	std::priority_queue<ROTNODE> queueRot;
	float v1, v2, v3, t, ct, ct2,st, st2;
	float tmp121, tmp122, tmp131, tmp132, tmp231, tmp232;
	float R11, R12, R13, R21, R22, R23, R31, R32, R33;
	std::map<int,std::vector<LINE3D>> inlier_from_rot;

	double remain_vol =  8 * PI * PI * PI;
	optInlier_ = -1;
	Eigen::Matrix3d R_init;
	optimizer_.optimize_R(RotationMSEThresh, lModel, pData, optR, R_init, optInlier_);
	int current_Inlier = 0;
	for(int i = 0; i < Nd; i++)
	{	
		int i_line_num = lModel[i].size();
		PLANE3D plane = pData[i];
		Eigen::Vector3d u = plane.u_;
		Eigen::Vector3d v = plane.v_;
		
		for(int j = 0; j < i_line_num; j++){
			Eigen::Vector3d d = lModel[i].at(j).d_;
			Eigen::Vector3d q1_R = (d.dot(R_init * u) * R_init * u + d.dot(R_init * v) * R_init * v).normalized();    
			double dot = d.dot(q1_R);
			
			double error;
			dot > 0 ? error = acos(dot) : error = M_PI - acos(dot);
			if(error < RotationMSEThresh) current_Inlier++;
		}
	}

	optR = R_init;

	queueRot.push(initNodeRot);
	long long count = 0;

	while(1)
	{	
		if(queueRot.empty())
		{
		  std::cout << "Rotation Queue Empty" << std::endl;
		  break;
		}

		// Access rotation cube with lowest lower bound...
		nodeRotParent = queueRot.top();
		// ...and remove it from the queue
		queueRot.pop();

		// Exit if the optError is less than or equal to the lower bound plus a small epsilon
		if(nodeRotParent.ub <= optInlier_)
		{	
			std::cout << nodeRotParent.ub << " " << optInlier_ << std::endl;
			break;
		}

		
		count ++;
		nodeRot.w = nodeRotParent.w/2;
		nodeRot.l = nodeRotParent.l+1;

		for(int j = 0; j < 8; j++)
		{	
			ub = 0;
			lb = 0;
			nodeRot.a = nodeRotParent.a + (j&1)*nodeRot.w;
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

			for(int i = 0; i<Nd; i++){                
				int line_num = lModel[i].size();
				Eigen::Vector3d u = pData.at(i).u_;
				Eigen::Vector3d v = pData.at(i).v_;
				for(int j = 0; j<line_num; j++){
					Eigen::Vector3d d = lModel[i].at(j).d_;
					Eigen::Vector3d q_i_R0_vec = (d.dot(nodeRot.SO3 *u) * nodeRot.SO3 * u + d.dot(nodeRot.SO3 *v) * nodeRot.SO3 * v).normalized(); 
					double dot = d.dot(q_i_R0_vec);
					double error;
					dot > 0 ? error = safe_acos(dot) : error = M_PI - safe_acos(dot);
					psi1[i].push_back(error);
					if(RotationMSEThresh - error > 0){
						lb++;
					}
				}
			}

			if(optInlier_ < 2*lb){
				optNodeRot = nodeRot;
				R_init = nodeRot.SO3;
				Eigen::Matrix3d R_optimized;
				int current_Inlier = optimizer_.optimize_R(RotationMSEThresh, lModel, pData, R_init, R_optimized, optInlier_);

				if(optInlier_ < current_Inlier)
				{	
					optInlier_ = current_Inlier;
					optR = R_optimized;

					std::map<int,std::vector<LINE3D>> current_inliers;
					for (int i=0; i< Nd; i++){
						int line_num = lModel[i].size();
						Eigen::Vector3d u = pData.at(i).u_;
						Eigen::Vector3d v = pData.at(i).v_;
						std::vector<LINE3D> temp_inliers;

						for(int j = 0 ; j<line_num;j++){
							Eigen::Vector3d d = lModel[i].at(j).d_;
							Eigen::Vector3d q_R = ((d.dot(optR * u)) * optR * u + (d.dot(optR * v)) * optR * v).normalized(); 
							double dot = d.dot(q_R);
							double error;
							dot > 0 ? error = acos(dot) : error = M_PI - acos(dot);
							if(error < RotationMSEThresh){
								temp_inliers.push_back(lModel[i].at(j));
							}
						}
						if(temp_inliers.size() > 0){
							current_inliers[i] = temp_inliers;
						}
					}
					inlier_from_rot = current_inliers;
				}

				std::priority_queue<ROTNODE> queueRotNew;

		        while(!queueRot.empty())
				{
					// RotNode node = queueRot.top();
					ROTNODE node = queueRot.top();
					queueRot.pop();
					if(node.ub > optInlier_)
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
				int line_num = lModel[i].size();
				for(int j = 0; j <line_num; j++){
					double d_R0_dist = psi1[i].at(j);
					if(RotationMSEThresh + current_maxRot - d_R0_dist > 0){
						ub++;
					} 
				}
			}

			for(int i = 0; i<Nd ; i++) psi1[i].clear();
			if(ub < optInlier_)
			{
				continue;
			}

    		nodeRot.ub = ub;
			nodeRot.lb = lb;
			queueRot.push(nodeRot);
		}
	}
	return inlier_from_rot;
}

void GoL2P::TranslationMCSearch(std::map<int,std::vector<LINE3D>> inlier_from_rot){
	
	float transX, transY, transZ;
    float lb, ub;
	TRANSNODE_MC nodeTrans, nodeTransParent;
    std::priority_queue<TRANSNODE_MC> queueTrans;

	std::vector<int> plane_ids;
    

	for (const auto& [key, value] : inlier_from_rot) {
        plane_ids.push_back(key);
    }


	std::map<int,std::vector<Eigen::Vector4d>> q2_fixed_values;
	int plane_num = inlier_from_rot.size();
	int inlier_num = 0;

	for(const auto & id : plane_ids){
		int line_num = inlier_from_rot[id].size();
		inlier_num += line_num;
		Eigen::Vector3d R_u = optR * pData.at(id).u_;
		Eigen::Vector3d R_v = optR* pData.at(id).v_;
		std::vector<Eigen::Vector4d> temp_q2;
		for(int j = 0; j<line_num; j++){
			Eigen::Vector3d d = inlier_from_rot[id].at(j).d_;
			Eigen::Vector3d q2 = (d.dot(R_v) * R_u - d.dot(R_u) * R_v).normalized();
			Eigen::Vector4d b_tilde = (inlier_from_rot[id].at(j).b_).homogeneous().normalized();
			Eigen::Vector4d q2_bar = Eigen::Vector4d::Zero();
						
			q2_bar.block(0,0,3,1) = q2;
			temp_q2.push_back(b_tilde.dot(q2_bar) * q2_bar);
		}
		q2_fixed_values[id] = temp_q2;
	}

	initNodeTrans_MC.ub = inlier_num;
	queueTrans.push(initNodeTrans_MC);

	long long count = 0;
	optInlier_ = -1;
	Eigen::Vector3d t_init;

	optimizer_.optimize_t_MC(TranslationMSEThresh,inlier_from_rot, pData, q2_fixed_values, plane_ids, optT, t_init, optR, optInlier_);
	int init_inlier = 0;


	for (int i = 0; i < plane_num; i++){
			int id = plane_ids.at(i);
			int line_num = inlier_from_rot[id].size();
			Eigen::Vector3d n = pData.at(id).normal_;
			Eigen::Vector3d x = pData.at(id).x1_;
			Eigen::Vector4d c_new_tilde= (n.dot(x-t_init) * optR * n).homogeneous().normalized();

			for (int j=0; j< line_num; j++){
				Eigen::Vector4d b_tilde = inlier_from_rot[id].at(j).b_.homogeneous().normalized();      
				Eigen::Vector4d q2_new = (q2_fixed_values[id].at(j) + b_tilde.dot(c_new_tilde) * c_new_tilde).normalized();
				double dot = b_tilde.dot(q2_new);
				double error;
				dot > 0 ? error = acos(dot) : error = M_PI - acos(dot);
				if(error < TranslationMSEThresh) init_inlier++;

			}
	}

	optT = t_init;
	
	while(1)
	{
		if(queueTrans.empty()){

		  std::cout << "Translation Queue Empty" << std::endl;
			break;
		}

		nodeTransParent = queueTrans.top();
		queueTrans.pop();

		if(nodeTransParent.ub <= optInlier_)
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

            std::map<int,std::vector<double>> current_distance;
            nodeTrans.x = nodeTransParent.x + (j&1)*nodeTrans.w ;
			nodeTrans.y = nodeTransParent.y + (j>>1&1)*nodeTrans.w ;
			nodeTrans.z = nodeTransParent.z + (j>>2&1)*nodeTrans.w ;
            
            transX = nodeTrans.x + nodeTrans.w/2;
			transY = nodeTrans.y + nodeTrans.w/2;
			transZ = nodeTrans.z + nodeTrans.w/2;

            Eigen::Vector3d t0 = {transX, transY, transZ};
            std::map<int,std::vector<double>> maxTransDist;
            for(int i = 0; i< plane_num; i++)
			{	
				int id = plane_ids.at(i);
				Eigen::Vector3d c_i_R_t0_vec = pData.at(id).normal_.dot(pData.at(id).x1_ - t0) * optR * pData.at(id).normal_;
				Eigen::Vector4d c_i_R_t0_tilde = c_i_R_t0_vec.homogeneous().normalized();
				int line_num = inlier_from_rot[id].size();
				std::vector<double> temp_distances;
				std::vector<Eigen::Vector4d> c_vert_values;

				 for(int k =0; k < 8; k++){
                    Eigen::Vector3d vertice;
                    vertice << nodeTrans.x + (k&1)*nodeTrans.w, nodeTrans.y + (k>>1&1)*nodeTrans.w, nodeTrans.z + (k>>2&1)*nodeTrans.w;
                    Eigen::Vector3d c_i_R_vert = pData.at(id).normal_.dot(pData.at(id).x1_ - vertice) * optR * pData.at(id).normal_;
					c_vert_values.push_back(c_i_R_vert.homogeneous().normalized());
                }
		
				std::vector<double> temp_psi_t;
				for(int j = 0; j<line_num; j++){
					Eigen::Vector4d b_tilde = inlier_from_rot[id].at(j).b_.homogeneous().normalized();
					Eigen::Vector4d q2_R_t0 = (q2_fixed_values[id].at(j) + b_tilde.dot(c_i_R_t0_tilde) * c_i_R_t0_tilde).normalized();
					temp_distances.push_back(ComputeGrassDist(q2_R_t0, b_tilde));
					double max_trans_dist = 0;
					for(const auto & c_vert : c_vert_values){
						Eigen::Vector4d q2_R_vert = (q2_fixed_values[id].at(j) + b_tilde.dot(c_vert) * c_vert).normalized();
						double dist = ComputeGrassDist(q2_R_t0,q2_R_vert);
						if(dist > max_trans_dist) max_trans_dist = dist;
					}
					if(max_trans_dist > M_PI /2) max_trans_dist = M_PI / 2;
					temp_psi_t.push_back(max_trans_dist);	  
				}

				current_distance[id] = temp_distances;
				maxTransDist[id]= temp_psi_t;

			}

            for(int i = 0; i < plane_num; i++)
            {	
				int id = plane_ids.at(i);
				int line_num = inlier_from_rot[id].size();
				for(int j = 0; j<line_num;j++){
					if(TranslationMSEThresh - current_distance[id].at(j) > 0){
						lb++;
					}
				}
            }
			


			if(optInlier_ < 2*lb)
			{
				// Update optimal error and rotation/translation nodes
				optNodeTrans_MC = nodeTrans;
				t_init << optNodeTrans_MC.x+optNodeTrans_MC.w/2, optNodeTrans_MC.y+optNodeTrans_MC.w/2, optNodeTrans_MC.z+optNodeTrans_MC.w/2;
			
				
				
				// Run ICP
				Eigen::Vector3d t_optimized;
				int current_Inlier = optimizer_.optimize_t_MC(TranslationMSEThresh, inlier_from_rot, pData, q2_fixed_values, plane_ids, t_init, t_optimized, optR, optInlier_);


				if(optInlier_ < current_Inlier)
				{	
					optInlier_ = current_Inlier;
					optT = t_optimized;
				}

				// Discard all rotation nodes with high lower bounds in the queue
				std::priority_queue<TRANSNODE_MC> queueTransNew;

		        while(!queueTrans.empty())
					{
					TRANSNODE_MC node = queueTrans.top();
					queueTrans.pop();
					if(node.ub > optInlier_)
						queueTransNew.push(node);
					else
						break;
				}

				queueTrans = queueTransNew;
			}


            for(int i = 0; i < plane_num; i++){
				int id = plane_ids.at(i);
				int line_num = inlier_from_rot[id].size();
				for(int j =0; j< line_num; j++){
	                double psi_temp = current_distance[id].at(j) - maxTransDist[id].at(j);
					if(TranslationMSEThresh - psi_temp > 0 ){
						ub++;
					}
				}
            }

			if(ub < optInlier_)
			{
				continue;
			}

    		nodeTrans.ub = ub;
			nodeTrans.lb = lb;
			queueTrans.push(nodeTrans);
		}
	if(count > 3000) break;

	}

}




double GoL2P::TranslationSearch(std::map<int,std::vector<LINE3D>> inlier_from_rot){
	
	std::vector<int> plane_ids;
    

	for (const auto& [key, value] : inlier_from_rot) {
        plane_ids.push_back(key);
    }

	std::map<int,std::vector<Eigen::Vector4d>> q2_fixed_values;
	int plane_num = inlier_from_rot.size();
	int inlier_num = 0;

	for(const auto & id : plane_ids){
		int line_num = inlier_from_rot[id].size();
		inlier_num += line_num;
		Eigen::Vector3d R_u = optR * pData.at(id).u_;
		Eigen::Vector3d R_v = optR* pData.at(id).v_;
		std::vector<Eigen::Vector4d> temp_q2;
		for(int j = 0; j<line_num; j++){
			Eigen::Vector3d d = inlier_from_rot[id].at(j).d_;
			Eigen::Vector3d q2 = (d.dot(R_v) * R_u - d.dot(R_u) * R_v).normalized();
			Eigen::Vector4d b_tilde = (inlier_from_rot[id].at(j).b_).homogeneous().normalized();
			Eigen::Vector4d q2_bar = Eigen::Vector4d::Zero();
			
			q2_bar.block(0,0,3,1) = q2;
			
			temp_q2.push_back(b_tilde.dot(q2_bar) * q2_bar);
		}
		q2_fixed_values[id] = temp_q2;
	}
	
	float transX, transY, transZ;
    float lb, ub, optErrorT, error;
	optErrorT = 1e3;
	TRANSNODE nodeTrans, nodeTransParent;
    std::priority_queue<TRANSNODE> queueTrans;
	queueTrans.push(initNodeTrans);
	long long count = 0;
	Eigen::Vector3d t_init;
	optErrorT = optimizer_.optimize_t(TranslationMSEThresh,inlier_from_rot, pData, q2_fixed_values, plane_ids, inlier_num, optT, t_init, optR, optErrorT);

	optT = t_init;
	TranslationSSEThresh = TranslationMSEThresh * inlier_num;

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
		// For inital division to 8 points
        for(int j = 0; j < 8; j++)
		{   
            std::map<int,std::vector<double>> current_distance;
            nodeTrans.x = nodeTransParent.x + (j&1)*nodeTrans.w ;
			nodeTrans.y = nodeTransParent.y + (j>>1&1)*nodeTrans.w ;
			nodeTrans.z = nodeTransParent.z + (j>>2&1)*nodeTrans.w ;
            
            transX = nodeTrans.x + nodeTrans.w/2;
			transY = nodeTrans.y + nodeTrans.w/2;
			transZ = nodeTrans.z + nodeTrans.w/2;

            Eigen::Vector3d t0 = {transX, transY, transZ};
            std::map<int,std::vector<double>> maxTransDist;
            for(int i = 0; i< plane_num; i++)
			{	
				int id = plane_ids.at(i);
				Eigen::Vector3d c_i_R_t0_vec = pData.at(id).normal_.dot(pData.at(id).x1_ - t0) * optR * pData.at(id).normal_;
				Eigen::Vector4d c_i_R_t0_tilde = c_i_R_t0_vec.homogeneous().normalized();
				int line_num = inlier_from_rot[id].size();
				std::vector<double> temp_distances;
				std::vector<Eigen::Vector4d> c_vert_values;
				 for(int k =0; k < 8; k++){
                    Eigen::Vector3d vertice;
                    vertice << nodeTrans.x + (k&1)*nodeTrans.w, nodeTrans.y + (k>>1&1)*nodeTrans.w, nodeTrans.z + (k>>2&1)*nodeTrans.w;
                    Eigen::Vector3d c_i_R_vert = pData.at(id).normal_.dot(pData.at(id).x1_ - vertice) * optR * pData.at(id).normal_;
					c_vert_values.push_back(c_i_R_vert.homogeneous().normalized());
                }
		
				std::vector<double> temp_psi_t;
				for(int j = 0; j<line_num; j++){
					Eigen::Vector4d b_tilde = inlier_from_rot[id].at(j).b_.homogeneous().normalized();
					Eigen::Vector4d q2_R_t0 = (q2_fixed_values[id].at(j) + b_tilde.dot(c_i_R_t0_tilde) * c_i_R_t0_tilde).normalized();
	                temp_distances.push_back(ComputeGrassDist(q2_R_t0, b_tilde));
					double max_trans_dist = 0;
					for(const auto & c_vert : c_vert_values){
						Eigen::Vector4d q2_R_vert = (q2_fixed_values[id].at(j) + b_tilde.dot(c_vert) * c_vert).normalized();
						double dist = ComputeGrassDist(q2_R_t0,q2_R_vert);
						if(dist > max_trans_dist) max_trans_dist = dist;
					}
					if(max_trans_dist > M_PI /2) max_trans_dist = M_PI / 2;

					temp_psi_t.push_back(max_trans_dist);	  
				}
				current_distance[id] = temp_distances;
				maxTransDist[id]= temp_psi_t;
			}

            ub = 0;
            for(int i = 0; i < plane_num; i++)
            {	
				int id = plane_ids.at(i);
				int line_num = inlier_from_rot[id].size();
				for(int j = 0; j<line_num;j++){

					ub += current_distance[id].at(j) * current_distance[id].at(j);
				}
            }

			if(ub < optErrorT)
			{
				// Update optimal error and rotation/translation nodes
				optNodeTrans = nodeTrans;

				// Run ICP
				t_init << optNodeTrans.x+optNodeTrans.w/2, optNodeTrans.y+optNodeTrans.w/2, optNodeTrans.z+optNodeTrans.w/2;
				Eigen::Vector3d t_optimized;
				error = optimizer_.optimize_t(TranslationMSEThresh,inlier_from_rot, pData, q2_fixed_values, plane_ids, inlier_num, t_init, t_optimized, optR, optErrorT);

				if(error < optErrorT)
				{	
					optErrorT = error;
					optT = t_optimized;
				}


				// Discard all rotation nodes with high lower bounds in the queue
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
            for(int i = 0; i < plane_num; i++){
				int id = plane_ids.at(i);
				int line_num = inlier_from_rot[id].size();
				for(int j =0; j< line_num; j++){
	                double psi_temp = current_distance[id].at(j) - maxTransDist[id].at(j);
					if(psi_temp < 0) psi_temp = 0;
	                lb += psi_temp * psi_temp;
				}
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
		if(count > 500) break;
	}
	return optErrorT;
}







std::pair<double,int> GoL2P::Register()
{	
    std::chrono::high_resolution_clock::time_point t_time_begin, t_time_end;
	std::chrono::microseconds t_diff_usec;  // For microseconds, use microseconds; for nanoseconds, use nanoseconds.
	Initialize();
    t_time_begin = std::chrono::high_resolution_clock::now();
	auto rotation_inliers = RotationSearch();
	
	double errorT = TranslationSearch(rotation_inliers);	
	// TranslationMCSearch(rotation_inliers);	

	t_time_end = std::chrono::high_resolution_clock::now();
	t_diff_usec = std::chrono::duration_cast<std::chrono::microseconds>(t_time_end - t_time_begin);
	double t_diff_sec = t_diff_usec.count() / 1000000.0;
	
	printf("time : %fms\n", (double)t_diff_sec * 1000);
	Clear();
	return std::make_pair(t_diff_sec *1000, optInlier_);
}
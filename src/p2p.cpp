#include <stdlib.h>
#include <math.h>
#include <chrono>
#include <stdio.h>
#include <vector>
#include "p2p.hpp"
#include "line_fitting.hpp"
#include "refinement.hpp"

using namespace pmc;

GoP2P::GoP2P()
{        
	initNodeRot.a = -PI;
	initNodeRot.b = -PI;
	initNodeRot.c = -PI;
	initNodeRot.w = 2*PI;
	initNodeRot.l = 0;
	initNodeRot.ub = Nd;

	initNodeRot_full.a = -PI;
	initNodeRot_full.b = -PI;
	initNodeRot_full.c = -PI;
	initNodeRot_full.w = 2*PI;
	initNodeRot_full.lb = 0;

	initNodeTrans.lb = 0;
}

void GoP2P::Initialize()
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


	// Initialise so-far-best rotation and translation nodes
	initNodeRot.ub = Nd;
	optNodeRot = initNodeRot;
	optNodeRot_full = initNodeRot_full;

	optNodeTrans = initNodeTrans;

	// Initialise so-far-best rotation and translation matrices
	optR = Eigen::Matrix3d::Identity();
	optT = Eigen::Vector3d::Zero();

	inlierNum = Nd;

	RotationSSEThresh = RotationMSEThresh * inlierNum;
}

void GoP2P::Clear()
{
	delete(psi1);

}

double GoP2P::InnerBnB(double * current_r_uncert, TRANSNODE* nodeTransOut, int level, Eigen::Matrix3d R0, double cur_maxRotAng){
	int i, j;
	float transX, transY, transZ;
	float lb, ub, optErrorT;
	float dis, maxTransDis;
	TRANSNODE nodeTrans, nodeTransParent;
	priority_queue<TRANSNODE> queueTrans;
	double * psi2 = new double[Nd];
	// Push top-level translation node into the priority queue
	queueTrans.push(initNodeTrans);
	optErrorT = 1e2;
	while(1){
		if(queueTrans.empty())
			break;

		nodeTransParent = queueTrans.top();
		queueTrans.pop();
		if(optErrorT-nodeTransParent.lb < RotationSSEThresh)
		{

			break;
		}

		nodeTrans.w = nodeTransParent.w/2;
		for(j = 0; j < 8; j++)
		{
			nodeTrans.x = nodeTransParent.x + (j&1)*nodeTrans.w ;
			nodeTrans.y = nodeTransParent.y + (j>>1&1)*nodeTrans.w ;
			nodeTrans.z = nodeTransParent.z + (j>>2&1)*nodeTrans.w ;

			transX = nodeTrans.x + nodeTrans.w/2;
			transY = nodeTrans.y + nodeTrans.w/2;
			transZ = nodeTrans.z + nodeTrans.w/2;
			Eigen::Vector3d t0;
			t0 << transX, transY, transZ;
			std::vector<double> maxTransDist;
			for(i = 0; i < Nd; i++)
			{

				Eigen::Vector4d tilde_c1_t0 = (pData.at(i).normal_.dot(pData.at(i).x1_ - t0) * pData.at(i).normal_).homogeneous().normalized();
				Eigen::Vector4d tilde_c2_r0 = (pModel.at(i).normal_.dot(pModel.at(i).x1_) * R0.transpose() * pModel.at(i).normal_).homogeneous().normalized();
				psi2[i] = ComputeGrassDist(tilde_c1_t0, tilde_c2_r0);

	
				
				if(current_r_uncert){

					psi1[i] -= cur_maxRotAng;
					psi2[i] -= current_r_uncert[i];
				}

				if(psi1[i] < 0){
					psi1[i] = 0;
				}  

                double max_trans_dist = 0;
                for(int k =0; k < 8; k++){
                    Eigen::Vector3d vertice;
                    vertice << nodeTrans.x + (k&1)*nodeTrans.w, nodeTrans.y + (k>>1&1)*nodeTrans.w, nodeTrans.z + (k>>2&1)*nodeTrans.w;
                    Eigen::Vector4d tilde_c_vert = (pData.at(i).normal_.dot(pData.at(i).x1_ - vertice) * pData.at(i).normal_).homogeneous().normalized();
					double dist = ComputeGrassDist(tilde_c1_t0,tilde_c_vert);
					if(dist > max_trans_dist) max_trans_dist = dist;
                }
                if(max_trans_dist > M_PI /2) max_trans_dist = M_PI / 2;
				maxTransDist.push_back(max_trans_dist);	     
			}


			ub = 0;
			for(i = 0; i < inlierNum; i++)
			{	
				double trans_part = psi2[i];
				if(trans_part < 0){
					trans_part = 0;
				}
				ub += psi1[i]* psi1[i] + trans_part * trans_part;
			}

			lb = 0;
			for(i = 0; i < inlierNum; i++)
			{
				double trans_part = psi2[i] - maxTransDist.at(i);
				if(trans_part < 0){
					trans_part = 0;
				}
				lb += psi1[i] * psi1[i] + trans_part * trans_part;
			}

			if(ub < optErrorT)
			{	
				optErrorT = ub;
				if(nodeTransOut)
					*nodeTransOut = nodeTrans;
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
	}
	delete psi2;
	return optErrorT;

}
double GoP2P::OuterBnB(){
	int i, j;
	ROTNODE_FULL nodeRot, nodeRotParent;
	TRANSNODE nodeTrans;
	float v1, v2, v3, t, ct, ct2,st, st2;
	float tmp121, tmp122, tmp131, tmp132, tmp231, tmp232;
	float R11, R12, R13, R21, R22, R23, R31, R32, R33;
	float lb, ub, error, dis;
	clock_t clockBeginICP;
	priority_queue<ROTNODE_FULL> queueRot;
	optError = 1e3;
	double * rot_uncert = new double[Nd];
	queueRot.push(initNodeRot_full);

	long long count = 0;
	while(1){
		if(queueRot.empty())
		{
		  cout << "Rotation Queue Empty" << endl;
		  cout << "Error*: " << optError << ", LB: " << lb << endl;
		  break;
		}

		// Access rotation cube with lowest lower bound...
		nodeRotParent = queueRot.top();
		// ...and remove it from the queue
		queueRot.pop();

		// Exit if the optError is less than or equal to the lower bound plus a small epsilon
		if((optError-nodeRotParent.lb) <= RotationSSEThresh)
		{
			cout << "Error*: " << optError << ", LB: " << nodeRotParent.lb << ", epsilon: " << RotationSSEThresh << endl;
			break;
		}

		if(count>0 && count%300 == 0)
			printf("LB=%f  L=%d\n",nodeRotParent.lb,nodeRotParent.l);
		count ++;
		
		// Subdivide rotation cube into octant subcubes and calculate upper and lower bounds for each
		nodeRot.w = nodeRotParent.w/2;
		nodeRot.l = nodeRotParent.l+1;
		for(j = 0; j < 8; j++){
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

			for(i = 0; i < Nd; i++)
			{
				Eigen::Vector3d R0_n1 = nodeRot.SO3 * pData.at(i).normal_;
				double dot = R0_n1.dot(pModel.at(i).normal_);
				double dist;
				dot > 0 ? dist = safe_acos(dot) : dist = PI - safe_acos(dot);
				psi1[i] = dist;
				rot_uncert[i] = safe_acos((pow(pData.at(i).normal_.dot(pData.at(i).x1_),2) * cos(maxRotAng[nodeRot.l]) +1) / (pow(pData.at(i).normal_.dot(pData.at(i).x1_),2) + 1));
			}

			double current_maxRot;
			if(maxRotAng[nodeRot.l] >= M_PI/2) current_maxRot = M_PI /2;
			else current_maxRot = maxRotAng[nodeRot.l];

			ub = InnerBnB(NULL /*Rotation Uncertainty Radius*/, &nodeTrans, nodeRot.l, nodeRot.SO3, current_maxRot);

			if(ub < optError){ 
				optError = ub;
				optNodeRot_full = nodeRot;
				optNodeTrans = nodeTrans;

				optR = nodeRot.SO3;
				optT << optNodeTrans.x+optNodeTrans.w/2, optNodeTrans.y+optNodeTrans.w/2, optNodeTrans.z+optNodeTrans.w/2;

				Eigen::Matrix3d R_optimized;
				Eigen::Vector3d t_optimized;

				error = optimizer_.optimize_R_t(pData, pModel,optR,R_optimized, optT, t_optimized, optError);
				if(error < optError)
				{
					optError = error;
					optR = R_optimized;
					optT = t_optimized;
					
					cout << "Error*: " << error << "(ICP " << (double)(clock() - clockBeginICP)/CLOCKS_PER_SEC << "s)" << endl;
				}

				priority_queue<ROTNODE_FULL> queueRotNew;
				while(!queueRot.empty())
				{
					ROTNODE_FULL node = queueRot.top();
					queueRot.pop();
					if(node.lb < optError)
						queueRotNew.push(node);
					else
						break;
				}
				queueRot = queueRotNew;
			}

			lb = InnerBnB(rot_uncert, NULL /*Translation Node*/, nodeRot.l, nodeRot.SO3, current_maxRot);
			if(lb >= optError)
			{
				continue;
			}

			// Update node and put it in queue
			nodeRot.ub = ub;
			nodeRot.lb = lb;
			queueRot.push(nodeRot);

		}
	}
	delete rot_uncert;
	return optError;
}

double GoP2P::RotationSearch_full(){
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

	optErrorR = optimizer_.optimize_R_full(RotationMSEThresh, pData, pModel, optR, R_init, optErrorR, Nm);
	

	optR = R_init;
	RotationSSEThresh = RotationMSEThresh * Nm;

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

		
	
		if(optErrorR-nodeRotParent.lb < RotationSSEThresh)
		{	
			std::cout << "OptErrorR: " << optErrorR << " nodeRotParent.lb " << nodeRotParent.lb << " SSE " << RotationSSEThresh << std::endl;
			break;
		}

		
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
			for(int i = 0; i<Nd ; i++){                
				Eigen::Vector3d n_i_R0_vec = nodeRot.SO3 * pData.at(i).normal_;
				Eigen::Vector3d pModel_n = pModel.at(i).normal_;
				double dot = n_i_R0_vec.dot(pModel_n);
				double error;
				dot > 0 ? error = safe_acos(dot) : error = PI - safe_acos(dot);
				psi1[i] = error;
				ub += error*error;

			}


			if(ub < optErrorR){
				optNodeRot_full = nodeRot;
				R_init = nodeRot.SO3;
				Eigen::Matrix3d R_optimized;
				error = optimizer_.optimize_R_full(RotationMSEThresh, pData, pModel, R_init, R_optimized, optErrorR, Nm);

				if(error < optErrorR)
				{	
					optErrorR = error;
					optR = R_optimized;

				}

				std::priority_queue<ROTNODE_FULL> queueRotNew;
		        while(!queueRot.empty())
				{
					ROTNODE_FULL node = queueRot.top();
					queueRot.pop();
					if(node.lb < optErrorR)
						queueRotNew.push(node);
					else
						break;
				}
				std::cout << queueRotNew.size() << std::endl;
				queueRot = queueRotNew;
			}
			lb = 0;

			for(int i = 0; i<Nd ; i++){       
				double current_maxRot;
				if(maxRotAng[nodeRot.l] >= M_PI/2) current_maxRot = M_PI /2;
				else current_maxRot = maxRotAng[nodeRot.l];          
				double psi_temp = psi1[i] - current_maxRot;
				if(psi_temp <0) psi_temp =0;
				lb += psi_temp * psi_temp;
			}


			if(lb >= optErrorR)
			{
				continue;
			}

    		nodeRot.ub = ub;
			nodeRot.lb = lb;
			queueRot.push(nodeRot);
		}
		count++;

		if(count > 1000) break;
	}

	return optErrorR;
}






std::vector<int> GoP2P::RotationSearch(){
	ROTNODE nodeRot, nodeRotParent;
	float lb, ub, error;
	clock_t clockBeginICL;
	
	std::priority_queue<ROTNODE> queueRot;
	float v1, v2, v3, t, ct, ct2,st, st2;
	float tmp121, tmp122, tmp131, tmp132, tmp231, tmp232;
	float R11, R12, R13, R21, R22, R23, R31, R32, R33;

	optInlier = -1;

    std::vector<int> inlier_from_rot;
	Eigen::Matrix3d R_init;
	optimizer_.optimize_R(pData, pModel, optR, R_init, RotationMSEThresh, optInlier);
	
	int current_Inlier = 0;

	optR = Eigen::Matrix3d::Identity();
	for(int i = 0; i < Nd; i++)
	{	
		std::vector<double> data_direction;
		Eigen::Vector3d pData_R = optR * pData[i].normal_;
		Eigen::Vector3d pModel_n = pModel[i].normal_;
		double dot = pData_R.dot(pModel_n);
		double error;
		dot > 0 ? error = safe_acos(dot) : error = PI - safe_acos(dot);
        if(error < RotationMSEThresh) current_Inlier++;
	}

	optInlier = current_Inlier;
	
	queueRot.push(initNodeRot);
	long long count = 0;

	while(1)
	{	
		if(queueRot.empty())
		{
		  break;
		}

		// Access rotation cube with lowest lower bound...
		nodeRotParent = queueRot.top();
		// ...and remove it from the queue
		queueRot.pop();

		if(nodeRotParent.ub <= optInlier)
		{	
			std::cout << nodeRotParent.ub << " " << optInlier << std::endl;

			break;
		}


		count ++;
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
				Eigen::Vector3d n_i_R0_vec = nodeRot.SO3 * pData.at(i).normal_;
				Eigen::Vector3d pModel_n = pModel.at(i).normal_;
				double dot = n_i_R0_vec.dot(pModel_n);
				double error;
				dot > 0 ? error = safe_acos(dot) : error = PI - safe_acos(dot);

				psi1[i] = error;

				if(RotationMSEThresh - psi1[i] > 0){
					lb++;
				}
			}

			if(optInlier < 2*lb){
				optNodeRot = nodeRot;
				R_init = nodeRot.SO3;
				Eigen::Matrix3d R_optimized;
				int current_Inlier = optimizer_.optimize_R(pData, pModel, R_init, R_optimized, RotationMSEThresh, optInlier);

				if(optInlier < current_Inlier)
				{	
					optInlier = current_Inlier;
					optR = R_optimized;
					// Update Correspondence Pair for MC 
                    std::vector<int> currnet_inlier_idx;
                    for(int i =0; i<Nd; i++){
                        Eigen::Vector3d n_i_R0_vec = optR * pData.at(i).normal_;
                        Eigen::Vector3d pModel_n = pModel.at(i).normal_;
                        double dot = n_i_R0_vec.dot(pModel_n);
                        double error;
                        dot > 0 ? error = safe_acos(dot) : error = PI - safe_acos(dot);
                        if(error < RotationMSEThresh) currnet_inlier_idx.push_back(i);
                    }
                    inlier_from_rot = currnet_inlier_idx;
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
				// std::cout << psi1[i] << " " <<  current_maxRot << std::endl;
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


        // Update Correspondence Pair for MC 


	return inlier_from_rot;
}



void GoP2P::TranslationSearch(std::vector<int> inlier_from_rot){
	std::vector<PLANE3D> inlier_data, inlier_model;
	int c_num = inlier_from_rot.size();

	for(const auto & idx : inlier_from_rot){
		inlier_data.push_back(pData.at(idx)); 
		inlier_model.push_back(pModel.at(idx));
	}

	float transX, transY, transZ;
    float lb, ub, optErrorT, error;
	TRANSNODE nodeTrans, nodeTransParent;
    std::priority_queue<TRANSNODE> queueTrans;
	queueTrans.push(initNodeTrans);
	long long count = 0;
	Eigen::Vector3d t_init;

	double updated_error = 0;

	optErrorT = optimizer_.optimize_t(inlier_data, inlier_model, optT, t_init, optR, optErrorT);
	optT = t_init;
	TranslationSSEThresh = TranslationMSEThresh * inlier_data.size();
	optErrorT = 1e3;
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
            std::vector<double> current_distance;
            nodeTrans.x = nodeTransParent.x + (j&1)*nodeTrans.w ;
			nodeTrans.y = nodeTransParent.y + (j>>1&1)*nodeTrans.w ;
			nodeTrans.z = nodeTransParent.z + (j>>2&1)*nodeTrans.w ;
            
            transX = nodeTrans.x + nodeTrans.w/2;
			transY = nodeTrans.y + nodeTrans.w/2;
			transZ = nodeTrans.z + nodeTrans.w/2;

            Eigen::Vector3d t0 = {transX, transY, transZ};
            std::vector<double> maxTransDist;


            for(int i = 0; i < c_num; i++)
			{	
				
				Eigen::Vector3d Ru = optR * inlier_data.at(i).u_;
				Eigen::Vector3d Rv = optR * inlier_data.at(i).v_;
				Eigen::Matrix<double,3,2> A; 
				A.block(0,0,3,1) = inlier_data.at(i).u_;
				A.block(0,1,3,1) = inlier_data.at(i).v_;

				Eigen::Vector4d Ru_bar, Rv_bar;
				Ru_bar << Ru(0), Ru(1), Ru(2), 0;
				Rv_bar << Rv(0), Rv(1), Rv(2), 0;

                Eigen::Vector4d d_i_R_t0_tilde = (inlier_data.at(i).normal_.dot(inlier_data.at(i).x1_ - t0) * optR * inlier_data.at(i).normal_).homogeneous().normalized();
				Eigen::Vector4d d_m_tilde = inlier_model.at(i).disp_.homogeneous().normalized();
				Eigen::Vector4d P_dm_t0 = ((Ru_bar.dot(d_m_tilde)) * Ru_bar + (Rv_bar.dot(d_m_tilde)) * Rv_bar + (d_i_R_t0_tilde.dot(d_m_tilde)) * d_i_R_t0_tilde).normalized();
				double dist = (P_dm_t0 - d_m_tilde).norm(); 
                current_distance.push_back(dist);
				
                double max_trans_dist = 0;
                for(int k =0; k < 8; k++){
                    Eigen::Vector3d vertice;
                    vertice << nodeTrans.x + (k&1)*nodeTrans.w, nodeTrans.y + (k>>1&1)*nodeTrans.w, nodeTrans.z + (k>>2&1)*nodeTrans.w;

					Eigen::Vector4d d_i_R_vertice_tilde = (optR * (Eigen::Matrix3d::Identity() - A * A.transpose()) * (inlier_data.at(i).x1_ - vertice)).homogeneous().normalized();

					Eigen::Vector4d P_dm_vertice = ((Ru_bar.dot(d_m_tilde)) * Ru_bar + (Rv_bar.dot(d_m_tilde)) * Rv_bar + (d_i_R_vertice_tilde.dot(d_m_tilde)) * d_i_R_vertice_tilde).normalized();

					double dist = (P_dm_t0 - P_dm_vertice).norm();
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
				error = optimizer_.optimize_t(inlier_data,inlier_model,t_init,t_optimized,optR, optErrorT);
				
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
		// std::cout << "Count: " << count << std::endl;
	}
}


double GoP2P::ComputeDispDist(Eigen::Vector3d disp1, Eigen::Vector3d disp2){
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

std::pair<double,int> GoP2P::Register()
{	

	Initialize();
    std::chrono::high_resolution_clock::time_point t_time_begin, t_time_end;
    std::chrono::microseconds t_diff_usec;
    t_time_begin = std::chrono::high_resolution_clock::now();


	auto rotation_inliers = RotationSearch();
	TranslationSearch(rotation_inliers);	

	// RotationSearch_full();
	// std::vector<int> full_id;
	// for(int i =0; i< Nd; i++){
	// 	full_id.push_back(i);
	// }
	// TranslationSearch(full_id);

	Clear();
    t_time_end = std::chrono::high_resolution_clock::now();
    t_diff_usec = std::chrono::duration_cast<std::chrono::microseconds>(t_time_end - t_time_begin);             
    double t_diff_sec = t_diff_usec.count() / 1000000.0;

	printf("time : %.6fs\n",t_diff_sec);
	Clear();

	return std::make_pair(t_diff_sec * 1000, optInlier);
}

/********************************************************************
Main Function for point cloud registration with Go-ICP Algorithm
Last modified: Feb 13, 2014

"Go-ICP: Solving 3D Registration Efficiently and Globally Optimally"
Jiaolong Yang, Hongdong Li, Yunde Jia
International Conference on Computer Vision (ICCV), 2013

Copyright (C) 2013 Jiaolong Yang (BIT and ANU)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

#include <iostream>
#include <fstream>
using namespace std;

#include "l2p.hpp"
#include "ConfigMap.hpp"

#define DEFAULT_OUTPUT_FNAME "output.txt"
#define DEFAULT_CONFIG_FNAME "config.txt"
#define DEFAULT_MODEL_FNAME "model.txt"
#define DEFAULT_DATA_FNAME "data.txt"

void parseInput(int argc, char **argv, string & modelFName, string & dataFName, string & configFName, string & outputFName);
void readConfig(string FName, GoL2P & gol2p);
void loadLines(string FName, int & N, std::map<int, std::vector<LINE3D>> & lines_);
void loadPlanes(string FName, int & N, std::vector<PLANE3D> & planes_);

int main(int argc, char** argv)
{
	int Nm, Nd;
	string modelFName, dataFName, configFName, outputFname;

	GoL2P gol2p;

	parseInput(argc, argv, modelFName, dataFName, configFName, outputFname);
	readConfig(configFName, gol2p);

	loadLines(modelFName, Nm, gol2p.lModel);
	loadPlanes(dataFName, Nd, gol2p.pData);
	
	gol2p.Nm = Nm;
	gol2p.Nd = Nd;

	cout << "Model ID: " << modelFName << " (" << gol2p.Nm << "), Data ID: " << dataFName << " (" << gol2p.Nd << ")" << endl;
	cout << "Registering..." << endl;
	auto [time, optInlierNum] = gol2p.Register();

	// cout << "Optimal Rotation Matrix:" << endl;
	// cout << gol2p.optR << endl;
	// cout << "Optimal Translation Vector:" << endl;
	// cout << gol2p.optT << endl;
	// cout << "Optimal Inlier Num:" << endl;
	// cout << optInlierNum << endl;

	ofstream ofile;
	ofile.open(outputFname.c_str(), std::ios_base::app);
	Eigen::Vector3d t = -gol2p.optR * gol2p.optT;
	Eigen::Quaterniond q(gol2p.optR);  // Convert a rotation matrix to a quaternion
	// ofile << t(0,0) << " " << t(1,0) << " " << t(2,0) << " "<< q.x() << " " << q.y() << " " << q.z() << " " << q.w() <<  "\n";

	ofile << gol2p.optR(0,0) << " " << gol2p.optR(0,1) << " " <<gol2p.optR(0,2) << " "<<
	gol2p.optR(1,0) << " " << gol2p.optR(1,1) << " " << gol2p.optR(1,2) << " " << gol2p.optR(2,0) << 
	" " << gol2p.optR(2,1) << " " << gol2p.optR(2,2) << " " << gol2p.optT(0,0) << " " << gol2p.optT(1,0) << 
	" " << gol2p.optT(2,0) << " " << time <<  "\n";
	ofile.close();
	return 0;
}

void parseInput(int argc, char **argv, string & modelFName, string & dataFName, string & configFName, string & outputFName)
{
	// Set default values
	modelFName = DEFAULT_MODEL_FNAME;
	dataFName = DEFAULT_DATA_FNAME;
	configFName = DEFAULT_CONFIG_FNAME;
	outputFName = DEFAULT_OUTPUT_FNAME;

	//cout << endl;
	//cout << "USAGE:" << "./GOICP <MODEL FILENAME> <DATA FILENAME> <NUM DOWNSAMPLED DATA POINTS> <CONFIG FILENAME> <OUTPUT FILENAME>" << endl;
	//cout << endl;
	if(argc != 5) return;
	modelFName = argv[1];
	dataFName = argv[2];
	configFName = argv[3];
	outputFName = argv[4];

	cout << "INPUT:" << endl;
	cout << "(modelFName)->(" << modelFName << ")" << endl;
	cout << "(dataFName)->(" << dataFName << ")" << endl;
	cout << "(configFName)->(" << configFName << ")" << endl;
	cout << "(outputFName)->(" << outputFName << ")" << endl;
	cout << endl;
}

void readConfig(string FName, GoL2P & gol2p)
{
	// Open and parse the associated config file
	ConfigMap config(FName.c_str());
	gol2p.initNodeRot.a = config.getF("rotMinX");
	gol2p.initNodeRot.b = config.getF("rotMinY");
	gol2p.initNodeRot.c = config.getF("rotMinZ");
	gol2p.initNodeRot.w = config.getF("rotWidth");

	gol2p.initNodeRot_full.a = config.getF("rotMinX");
	gol2p.initNodeRot_full.b = config.getF("rotMinY");
	gol2p.initNodeRot_full.c = config.getF("rotMinZ");
	gol2p.initNodeRot_full.w = config.getF("rotWidth");


	gol2p.RotationMSEThresh = config.getF("RotationMSEThresh");
	gol2p.TranslationMSEThresh = config.getF("TranslationMSEThresh");

	gol2p.initNodeTrans.x = config.getF("transMinX");
	gol2p.initNodeTrans.y = config.getF("transMinY");
	gol2p.initNodeTrans.z = config.getF("transMinZ");
	gol2p.initNodeTrans.w = config.getF("transWidth");


	gol2p.initNodeTrans_MC.x = config.getF("transMinX");
	gol2p.initNodeTrans_MC.y = config.getF("transMinY");
	gol2p.initNodeTrans_MC.z = config.getF("transMinZ");
	gol2p.initNodeTrans_MC.w = config.getF("transWidth");



	cout << "CONFIG:" << endl;
	config.print();
	cout << endl;
}

void loadLines(string FName, int & N, std::map<int, std::vector<LINE3D>> & lines_){
	int i;
	ifstream ifile;

	ifile.open(FName.c_str(), ifstream::in);
	if(!ifile.is_open())
	{
		cout << "Unable to open point file '" << FName << "'" << endl;
		exit(-1);
	}
	ifile >> N; // First line has number of points to follow
    std::vector<LINE3D> temp_lines;
    int id = 0;
	for(i = 0; i < N; i++)
	{
		double sp_x, sp_y, sp_z, ep_x, ep_y, ep_z;
        int id_new = 0;
		ifile >> sp_x >> sp_y >> sp_z >> ep_x >> ep_y >> ep_z >> id_new;

        if(id_new > id){
            lines_[id] = temp_lines;
            id++;
            temp_lines.clear();
        }

		Eigen::Vector3d sp, ep;
		sp << sp_x, sp_y, sp_z;
		ep << ep_x, ep_y, ep_z;

		LINE3D line(sp,ep,id_new);
		
		// std::cout << "Line d : " << line.d_ << std::endl;
		// std::cout << "Line b : " << line.b_ << std::endl;

		temp_lines.push_back(line);

        if(i == N-1){
            lines_[id_new] = temp_lines;
            temp_lines.clear();
        }
	}

	ifile.close();

}


void loadPlanes(string FName, int & N, std::vector<PLANE3D> & planes_){
	int i;
	ifstream ifile;

	ifile.open(FName.c_str(), ifstream::in);
	if(!ifile.is_open())
	{
		cout << "Unable to open point file '" << FName << "'" << endl;
		exit(-1);
	}
	ifile >> N; // First line has number of points to follow
	for(i = 0; i < N; i++)
	{
		double p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z;
        int id;
		ifile >> p1_x >> p1_y >> p1_z >> p2_x >> p2_y >> p2_z >> p3_x >> p3_y >> p3_z >> id;

		Eigen::Vector3d x1, x2, x3;
		x1 << p1_x, p1_y, p1_z;
		x2 << p2_x, p2_y, p2_z;
        x3 << p3_x, p3_y, p3_z;
		PLANE3D plane(x1,x2,x3,id);
		// std::cout << "Plane u: " << plane.u_ << std::endl;
		// std::cout << "Plane v: " << plane.v_ << std::endl;

		planes_.push_back(plane);
	}

	ifile.close();

}
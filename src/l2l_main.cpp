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

#include "l2l.hpp"
#include "ConfigMap.hpp"


void parseInput(int argc, char **argv, string & modelFName, string & dataFName, string & configFName, string & outputFName);
void readConfig(string FName, GoL2L & goL2L);
void loadLines(string FName, int & N,std::vector<LINE3D> & lines_, std::vector<Eigen::Vector3d> & dir_, std::vector<Eigen::Vector3d> & disp_);

int main(int argc, char** argv)
{
	int Nm, Nd;
	string modelFName, dataFName, configFName, outputFname;
	std::vector<LINE3D> lModel, lData;
	GoL2L goL2L;

	parseInput(argc, argv, modelFName, dataFName, configFName, outputFname);
	readConfig(configFName, goL2L);
	std::cout << configFName << std::endl;

	
	loadLines(modelFName, Nm, goL2L.lModel, goL2L.d_Model, goL2L.b_Model);
	loadLines(dataFName, Nd, goL2L.lData, goL2L.d_Data, goL2L.b_Data);
	
	goL2L.Nm = Nm;
	goL2L.Nd = Nd;


	cout << "Model ID: " << modelFName << " (" << goL2L.Nm << "), Data ID: " << dataFName << " (" << goL2L.Nd << ")" << endl;
	cout << "Registering..." << endl;
	auto [time, opterror] = goL2L.Register();

	cout << "Optimal Rotation Matrix:" << endl;
	cout << goL2L.optR << endl;
	cout << "Optimal Translation Vector:" << endl;
	cout << goL2L.optT << endl;

	ofstream ofile;
	ofile.open(outputFname.c_str(), std::ios_base::app);
	Eigen::Vector3d t = -goL2L.optR * goL2L.optT;
	Eigen::Quaterniond q(goL2L.optR);  // Convert a rotation matrix to a quaternion
	ofile << t(0,0) << " " << t(1,0) << " " << t(2,0) << " "<< q.x() << " " << q.y() << " " << q.z() << " " << q.w() <<  "\n";


	ofile.close();

	return 0;
}

void parseInput(int argc, char **argv, string & modelFName, string & dataFName, string & configFName, string & outputFName)
{
	// Set default values

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

void readConfig(string FName, GoL2L & goL2L)
{
	// Open and parse the associated config file
	std::cout << FName << std::endl;
	ConfigMap config(FName.c_str());
	goL2L.initNodeRot.a = config.getF("rotMinX");
	goL2L.initNodeRot.b = config.getF("rotMinY");
	goL2L.initNodeRot.c = config.getF("rotMinZ");
	goL2L.initNodeRot.w = config.getF("rotWidth");

	goL2L.RotationMSEThresh = config.getF("RotationMSEThresh");
	goL2L.TranslationMSEThresh = config.getF("TranslationMSEThresh");

	goL2L.initNodeTrans.x = config.getF("transMinX");
	goL2L.initNodeTrans.y = config.getF("transMinY");
	goL2L.initNodeTrans.z = config.getF("transMinZ");
	goL2L.initNodeTrans.w = config.getF("transWidth");

	goL2L.initNodeTrans_MC.x = config.getF("transMinX");
	goL2L.initNodeTrans_MC.y = config.getF("transMinY");
	goL2L.initNodeTrans_MC.z = config.getF("transMinZ");
	goL2L.initNodeTrans_MC.w = config.getF("transWidth");

	// If < 0.1% trimming specified, do no trimming

	cout << "CONFIG:" << endl;
	config.print();
	//cout << "(doTrim)->(" << goL2L.doTrim << ")" << endl;
	cout << endl;
}

void loadLines(string FName, int & N,std::vector<LINE3D> & lines_, std::vector<Eigen::Vector3d> & dir_, std::vector<Eigen::Vector3d> & disp_){
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
		double sp_x, sp_y, sp_z, ep_x, ep_y, ep_z;
		ifile >> sp_x >> sp_y >> sp_z >> ep_x >> ep_y >> ep_z;
		Eigen::Vector3d sp, ep;
		sp << sp_x, sp_y, sp_z;
		ep << ep_x, ep_y, ep_z;
		LINE3D line(sp,ep);
		lines_.push_back(line);
		dir_.push_back(line.d_);
		disp_.push_back(line.b_);
	}

	ifile.close();

}
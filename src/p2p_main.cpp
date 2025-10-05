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

#include "p2p.hpp"
#include "ConfigMap.hpp"

#define DEFAULT_OUTPUT_FNAME "output.txt"
#define DEFAULT_CONFIG_FNAME "config.txt"
#define DEFAULT_MODEL_FNAME "model.txt"
#define DEFAULT_DATA_FNAME "data.txt"

void parseInput(int argc, char **argv, string & modelFName, string & dataFName, string & configFName, string & outputFName);
void readConfig(string FName, GoP2P & gop2p);
void loadPlanes(string FName, int & N, std::vector<PLANE3D> & planes_);

int main(int argc, char** argv)
{
    int Nm, Nd;
	string modelFName, dataFName, configFName, outputFname;
	std::vector<PLANE3D> lModel, lData;
	GoP2P gop2p;

	parseInput(argc, argv, modelFName, dataFName, configFName, outputFname);
	readConfig(configFName, gop2p);
	std::cout << configFName << std::endl;

	
	loadPlanes(modelFName, Nm, gop2p.pModel);
	loadPlanes(dataFName, Nd, gop2p.pData);
	
	gop2p.Nm = Nm;
	gop2p.Nd = Nd;


	cout << "Model ID: " << modelFName << " (" << gop2p.Nm << "), Data ID: " << dataFName << " (" << gop2p.Nd << ")" << endl;
	cout << "Registering..." << endl;
	auto [time, opterror] = gop2p.Register();



	ofstream ofile;
	ofile.open(outputFname.c_str(), std::ios_base::out);
	Eigen::Vector3d t = -gop2p.optR * gop2p.optT;
	Eigen::Quaterniond q(gop2p.optR);  // Convert a rotation matrix to a quaternion


	ofile << gop2p.optR(0,0) << " " << gop2p.optR(0,1) << " " <<gop2p.optR(0,2) << " "<<
	gop2p.optR(1,0) << " " << gop2p.optR(1,1) << " " << gop2p.optR(1,2) << " " << gop2p.optR(2,0) << 
	" " << gop2p.optR(2,1) << " " << gop2p.optR(2,2) << " " << gop2p.optT(0,0) << " " << gop2p.optT(1,0) << 
	" " << gop2p.optT(2,0) << " " << time <<  "\n";

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

void readConfig(string FName, GoP2P & gop2p)
{
	// Open and parse the associated config file
	ConfigMap config(FName.c_str());
	gop2p.initNodeRot.a = config.getF("rotMinX");
	gop2p.initNodeRot.b = config.getF("rotMinY");
	gop2p.initNodeRot.c = config.getF("rotMinZ");
	gop2p.initNodeRot.w = config.getF("rotWidth");

	gop2p.RotationMSEThresh = config.getF("RotationMSEThresh");
	gop2p.TranslationMSEThresh = config.getF("TranslationMSEThresh");

	gop2p.initNodeTrans.x = config.getF("transMinX");
	gop2p.initNodeTrans.y = config.getF("transMinY");
	gop2p.initNodeTrans.z = config.getF("transMinZ");
	gop2p.initNodeTrans.w = config.getF("transWidth");




	cout << "CONFIG:" << endl;
	config.print();
	cout << endl;
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


		planes_.push_back(plane);
	}

	ifile.close();

}
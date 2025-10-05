#include "refinement.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>

void loadLines(std::string FName, int & N, std::vector<std::vector<double>> & lines_){
	int i;
	std::ifstream ifile;

	ifile.open(FName.c_str(), std::ifstream::in);
	if(!ifile.is_open())
	{
		std::cout << "Unable to open point file '" << FName << "'" << std::endl;
		exit(-1);
	}
	ifile >> N; // First line has number of points to follow
	for(i = 0; i < N; i++)
	{
		double sp_x, sp_y, sp_z, ep_x, ep_y, ep_z;
        int id_new = 0;
		ifile >> sp_x >> sp_y >> sp_z >> ep_x >> ep_y >> ep_z >> id_new;

        std::vector<double> line = {sp_x, sp_y, sp_z, ep_x, ep_y, ep_z};

		lines_.push_back(line);

	}

	ifile.close();
}

void loadInlierLines(std::string FName, int & N, std::vector<std::vector<double>> & lines_, std::vector<int> inlier_idx){
	int i;
	std::ifstream ifile;

	ifile.open(FName.c_str(), std::ifstream::in);
	if(!ifile.is_open())
	{
		std::cout << "Unable to open point file '" << FName << "'" << std::endl;
		exit(-1);
	}
	ifile >> N; // First line has number of points to follow
	for(i = 0; i < N; i++)
	{
		if(inlier_idx[i] == 1){
			double sp_x, sp_y, sp_z, ep_x, ep_y, ep_z;
			int id_new = 0;
			ifile >> sp_x >> sp_y >> sp_z >> ep_x >> ep_y >> ep_z >> id_new;

			std::vector<double> line = {sp_x, sp_y, sp_z, ep_x, ep_y, ep_z};

			

			lines_.push_back(line);
		}

	}

	ifile.close();
}


void loadPlanes(std::string FName, int & N, std::vector<std::vector<double>> & planes_){
	int i;
	std::ifstream ifile;

	ifile.open(FName.c_str(), std::ifstream::in);
	if(!ifile.is_open())
	{
		std::cout << "Unable to open point file '" << FName << "'" << std::endl;
		exit(-1);
	}
	ifile >> N; // First line has number of points to follow
	for(i = 0; i < N; i++)
	{
		double p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z;
        int id;
		ifile >> p1_x >> p1_y >> p1_z >> p2_x >> p2_y >> p2_z >> p3_x >> p3_y >> p3_z >> id;

        std::vector<double> plane = {p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z};


		planes_.push_back(plane);
	}
	ifile.close();
}

void loadInlierPlanes(std::string FName, int & N, std::vector<std::vector<double>> & planes_, std::vector<int> inlier_idx){
	int i;
	std::ifstream ifile;

	ifile.open(FName.c_str(), std::ifstream::in);
	if(!ifile.is_open())
	{
		std::cout << "Unable to open point file '" << FName << "'" << std::endl;
		exit(-1);
	}
	ifile >> N; // First line has number of points to follow
	for(i = 0; i < N; i++)
	{
		if(inlier_idx[i] == 1){
			double p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z;
			int id;
			ifile >> p1_x >> p1_y >> p1_z >> p2_x >> p2_y >> p2_z >> p3_x >> p3_y >> p3_z >> id;

			std::vector<double> plane = {p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z};
			// std::cout << "Plane u: " << plane.u_ << std::endl;
			// std::cout << "Plane v: " << plane.v_ << std::endl;

			planes_.push_back(plane);
		}

	}
	ifile.close();
}


Eigen::Matrix3d to_SO3(const Eigen::Vector3d &v)
{
    double theta = v.norm();
    if (theta == double(0))
        return Eigen::Matrix3d::Identity();


    Eigen::Matrix3d v_hat;
	v_hat << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
    Eigen::Matrix3d mat = Eigen::Matrix3d::Identity() + (sin(theta) / theta) * v_hat + ((1.0 - cos(theta)) / (theta * theta)) * (v_hat * v_hat);

    return mat;
}

template <typename T>
bool loadMatrix(const std::string& filename, std::vector<std::vector<T>>& matrix) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open the file " << filename << std::endl;
        return false;
    }
    
    std::string line;
    while (getline(file, line)) {
        std::stringstream ss(line);
        std::vector<T> row;
        T value;
        
        while (ss >> value) {
            row.push_back(value);
        }
        
        matrix.push_back(row);
    }

    file.close();
    return true;
}


int main(int argc, char** argv){
	// checkerboard refine
	std::vector<std::string> planeFs;
	planeFs.push_back("/home/jaeho/random_pnl/chessboard/plane1.txt");
	planeFs.push_back("/home/jaeho/random_pnl/chessboard/plane2.txt");
	planeFs.push_back("/home/jaeho/random_pnl/chessboard/plane5.txt");
	planeFs.push_back("/home/jaeho/random_pnl/chessboard/plane6.txt");
	planeFs.push_back("/home/jaeho/random_pnl/chessboard/plane13.txt");
	planeFs.push_back("/home/jaeho/random_pnl/chessboard/plane14.txt");
	planeFs.push_back("/home/jaeho/random_pnl/chessboard/plane15.txt");
	planeFs.push_back("/home/jaeho/random_pnl/chessboard/plane17.txt");
	planeFs.push_back("/home/jaeho/random_pnl/chessboard/plane18.txt");
	planeFs.push_back("/home/jaeho/random_pnl/chessboard/plane19.txt");
    std::ifstream file("/home/jaeho/Downloads/est_pose_s_asp3l.txt");  // Open the file
	std::vector<Eigen::Matrix4d> poses;
	std::string data;
	if (file.is_open()) {
			// Read file line by line
			while (std::getline(file, data)) {
				std::stringstream ss(data);  // Use stringstream to split the line into values
				double value;
				std::vector<double> row;  // To store the current row

				// Read each number in the line and store it in the row vector
				while (ss >> value) {
					row.push_back(value);
				}
				Eigen::Vector3d rotvec = {row.at(3), row.at(4), row.at(5)};
				Eigen::Vector3d trans = {row.at(0), row.at(1), row.at(2)};
				Eigen::Matrix3d rotmat = to_SO3(rotvec);
				Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
				T.block(0,0,3,3) = rotmat.transpose();
				T.block(0,3,3,1) = -rotmat.transpose() *trans;
				poses.push_back(T);
				// Add the row to the data vector
			}
	}
			// Close the file
	file.close();
    int pair_num;
    std::vector<std::vector<double>> lines;
	std::string lineF = "/home/jaeho/random_pnl/chessboard/lines.txt";
	loadLines(lineF, pair_num, lines);
	double mean_time = 0;
	for(int i = 0; i <planeFs.size(); i++){
    	int plane_num;
		std::vector<std::vector<double>> planes;
		std::string planeF = planeFs.at(i);
		loadPlanes(planeF, plane_num, planes);
		Eigen::Matrix3d R_opt; 
		Eigen::Vector3d t_opt;

		
		std::chrono::high_resolution_clock::time_point t_time_begin, t_time_end;
		std::chrono::microseconds t_diff_usec;  // For microseconds, use microseconds; for nanoseconds, use nanoseconds.

		t_time_begin = std::chrono::high_resolution_clock::now();
		auto result = L2Poptimize(lines, planes, poses.at(i));
		t_time_end = std::chrono::high_resolution_clock::now();
		t_diff_usec = std::chrono::duration_cast<std::chrono::microseconds>(t_time_end - t_time_begin);
		
		double t_diff_msec = t_diff_usec.count() / 1000.0;
		mean_time += t_diff_msec;

		R_opt = result.first.transpose();
		double theta = acos((R_opt.trace() - 1.0) / 2.0);
		Eigen::Vector3d axis;
		axis << R_opt(2, 1) - R_opt(1, 2), R_opt(0, 2) - R_opt(2, 0), R_opt(1, 0) - R_opt(0, 1);
		axis.normalize();
		Eigen::Vector3d rotvec_opt = theta * axis;
		t_opt = -result.first.transpose() * result.second;
		std::cout << t_opt(0) << "\t" << t_opt(1) << "\t" << t_opt(2)<< "\t" << rotvec_opt(0) << "\t" << rotvec_opt(1) << "\t" << rotvec_opt(2) << std::endl;
	}
	std::cout << mean_time / 10 << std::endl;


	return 0;
}
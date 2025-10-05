# GrassmannRegistration (ICCV25, <em>highlight</em>)

This repository contains the code for GrassmannRegistration, an official implementation for the paper **Registration beyond Points: General Affine Subspace Alignment via Geodesic Distance on Grassmann Manifold**, which is accepted to ICCV 2025. GrassmannRegistration leverages geodesic distance on Grassmann manifold for cost function, leading to better convergence and extension to global optimization algorithm.

[[arXiv]](https://arxiv.org/abs/2507.17998) [[BibTex]](#bibtex)

GrassmannRegistration Supports
   
- Line to line registration
- Line to plane registration
- Plane to plane registration

## News
- 25.10.05 We released an inital version of GrassmannRegistration! 

## 1. Dependencies
The code was tested on following dependencies:
- Ubuntu 20.04
- Ceres Solver (2.2.0)
- Eigen3
- pmc (Parallel Maximum Clique Library)
## 2. Install

```bash
git clone https://github.com/joomeok/GrassmannRegistration.git
cd GrassmannRegistration
git submodule update --init --recursive
mkdir build
cd build 
cmake .. && make -j4
```
## 3. Examples
For now, we provide readily executable example only for BnB solver. However, you can refine your initial pose obtained from another solver with `refine_main.cpp`. Refer to detailed implementation within it for your specific example.

### 3.1 Data Format
For current version, you need two input files (features with respect to model and data frames, respectively). The format is:

### Line feature
```txt
$(Feature Num)
sp_x sp_y sp_z ep_x ep_y ep_z id
...
```
### Plane feature
```txt
$(Feature Num)
p1_x p1_y p1_z p2_x p2_y p2_z p3_x p3_y p3_z id
...
```
For line feature, `sp` and `ep` denote start point and end point respectively, while `p1`, `p2`, and `p3` denote 3 linearly independent points on plane feature. `id` contains the correspondence information. 

### 3.2 RGB-D Odometry
Run:
```bash
./RGBD_odom.sh
```
Refer to `toy_example.ipynb` for stacking delta pose to plot trajectory. You can compare resulted trajectory with ground truth trajectory in `/data` folder.

### 3.3 Line to Plane Registration

Refer to data generation code within `toy_example.ipynb`. You can run registration with:
```bash
./GoL2P $(YourDirectory)/lines.txt $(YourDirectory)/plane.txt ../config.txt $(YourDirectory)/result.txt
```

### 3.4 Custom Data
You can solve your own registration problem by running one of three binary files: `GoL2L`, `GoL2P`, and `GoP2P`. If you want to solve your plane-to-plane registration, try this command:

```bash
./GoP2P $(YourDirectory)/model_plane.txt $(YourDirectory)/data_plane.txt ../config.txt $(YourDirectory)/result.txt
```

## Acknowledgements
This project is based on code from [yangjialong/Go-ICP](https://github.com/yangjiaolong/Go-ICP).

## BibTex
```
@INPROCEEDINGS{jhshin-2025-iccv,  
    AUTHOR = { Jaeho Shin and Hyeonjae Gil and Junwoo Jang and Maani Ghaffari and Ayoung Kim },  
    TITLE = { Registration beyond Points: General Affine Subspace Alignment via Geodesic Distance on Grassmann Manifold },  
    BOOKTITLE = { Proceedings of the IEEE/CVF International Conference on Computer Vision (ICCV) },  
    YEAR = { 2025 },  
    MONTH = { October },  
    ADDRESS = { Honolulu, Hawaii }
}
```


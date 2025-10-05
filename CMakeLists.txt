cmake_minimum_required(VERSION 3.1)
project(GrassRegist)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC -fopenmp")

set(IO_SRC
    3rdparty/pmc/pmc_utils.cpp
    3rdparty/pmc/pmc_graph.cpp
    3rdparty/pmc/pmc_clique_utils.cpp
)
set(PMC_SRC
    3rdparty/pmc/pmc_heu.cpp
    3rdparty/pmc/pmc_maxclique.cpp
    3rdparty/pmc/pmcx_maxclique.cpp
    3rdparty/pmc/pmcx_maxclique_basic.cpp
)
set(BOUND_LIB_SRC
    3rdparty/pmc/pmc_cores.cpp
)

find_package(Ceres REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)

include_directories(
    ${CERES_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    include
    src
    3rdparty/pmc
)

add_library(libpmc SHARED
    ${IO_SRC}
    ${PMC_SRC}
    ${BOUND_LIB_SRC}
    3rdparty/pmc/pmc_lib.cpp
)

add_executable(GoICL
    src/l2l_main.cpp
    src/l2l.cpp
    src/line_fitting.cpp
    src/ConfigMap.cpp
    src/StringTokenizer.cpp
    src/line_3d.cpp
    src/util_funcs.cpp
)

add_executable(GoL2P
    src/l2p_main.cpp
    src/l2p.cpp
    src/l2p_fitting.cpp
    src/ConfigMap.cpp
    src/StringTokenizer.cpp
    src/line_3d.cpp
    src/plane_3d.cpp
    src/util_funcs.cpp
    src/refinement.cpp
)

add_executable(GoP2P
    src/p2p_main.cpp
    src/p2p.cpp
    src/p2p_fitting.cpp
    src/ConfigMap.cpp
    src/StringTokenizer.cpp
    src/plane_3d.cpp
    src/util_funcs.cpp
    src/refinement.cpp
)

add_executable(Refinement
    src/refine_main.cpp
    src/refinement.cpp
)

target_link_libraries(GoICL libpmc ${CERES_LIBRARIES} Eigen3::Eigen)
target_link_libraries(GoL2P ${CERES_LIBRARIES} Eigen3::Eigen)
target_link_libraries(GoP2P ${CERES_LIBRARIES} Eigen3::Eigen)
target_link_libraries(Refinement ${CERES_LIBRARIES} Eigen3::Eigen)


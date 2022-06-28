cmake_minimum_required(VERSION 3.15)
project(BI-GMRF)

set(CMAKE_CXX_STANDARD 14)

include_directories(/usr/local/include/eigen3 .)

add_executable(BI-GMRF main.cpp Constant.cpp Constant.h Generate.cpp Generate.h rnorm.cpp runif.cpp runifVec.cpp rgamma.cpp strToVec.cpp Initial.cpp Initial.h Mcmc.cpp Mcmc.h)

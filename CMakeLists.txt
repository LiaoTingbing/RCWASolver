

cmake_minimum_required(VERSION 3.10 )

# 设置工具链
set(CMAKE_TOOLCHAIN_FILE
        E:/VisualStudioFiles/vcpkg/vcpkg/vcpkg/scripts/buildsystems/vcpkg.cmake
)

project(rcwa )

add_executable(
        ${PROJECT_NAME} include/common.h include/convulationmatrix.h include/fftmatlab.h include/fftshiftmatlab.h include/rcwa.h include/userfunc.h
        src/convulationmatrix.cpp src/fftmatlab.cpp src/fftshiftmatlab.cpp src/rcwa.cpp src/userfunc.cpp
        examples1/example1.cpp
)


find_package(Armadillo CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE armadillo)


# this is heuristically generated, and may not be correct
find_package(FFTW3 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE FFTW3::fftw3)

find_package(FFTW3f CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE FFTW3::fftw3f)

find_package(FFTW3l CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE FFTW3::fftw3l)
 
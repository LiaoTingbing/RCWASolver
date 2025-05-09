

cmake_minimum_required(VERSION 3.10 )

set(CMAKE_TOOLCHAIN_FILE
E:/VisualStudioFiles/vcpkg/vcpkg/scripts/buildsystems/vcpkg.cmake)

project(rcwa )

add_executable(
${PROJECT_NAME}
	"include/common.h" "include/convulationmatrix.h" "include/fftmatlab.h" "include/fftshiftmatlab.h" "include/rcwa.h" "include/userfunc.h"
"src/convulationmatrix.cpp" "src/fftmatlab.cpp" "src/fftshiftmatlab.cpp" "examples1/example1.cpp" "src/rcwa.cpp" "src/userfunc.cpp"
)

 find_package(OpenMP REQUIRED)
 target_link_libraries(rcwa
  PUBLIC
      OpenMP::OpenMP_CXX
 )


target_include_directories(
	rcwa PRIVATE
	E:/VisualStudioFiles/vcpkg/ku/armadillo-14.4.1/include
)
target_link_libraries(
	rcwa PRIVATE
	E:/VisualStudioFiles/vcpkg/ku/armadillo-14.4.1/examples/lib_win64/libopenblas.lib
)

find_package(FFTW3 CONFIG REQUIRED)
target_link_libraries(    ${PROJECT_NAME}
PRIVATE FFTW3::fftw3)

find_package(FFTW3f CONFIG REQUIRED)
target_link_libraries(    ${PROJECT_NAME}
PRIVATE FFTW3::fftw3f)

find_package(FFTW3l CONFIG REQUIRED)
target_link_libraries(    ${PROJECT_NAME}
PRIVATE FFTW3::fftw3l)
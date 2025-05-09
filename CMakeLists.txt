

cmake_minimum_required(VERSION 3.10 )

project(rcwa )

add_executable(
${PROJECT_NAME}
	include/common.h "include/convulationmatrix.h" "include/fftmatlab.h" "include/fftshiftmatlab.h" "include/rcwa.h" "include/userfunc.h"
"src/convulationmatrix.cpp" "src/fftmatlab.cpp" "src/fftshiftmatlab.cpp" "examples1/example1.cpp" "src/rcwa.cpp" "src/userfunc.cpp"
)


#find_package(Armadillo CONFIG REQUIRED)
#target_link_libraries(${PROJECT_NAME} PRIVATE armadillo)

target_include_directories(
	${PROJECT_NAME} PRIVATE
	E:/VisualStudioFiles/vcpkg/ku/armadillo-14.4.1/include
)
target_link_libraries(
	${PROJECT_NAME}
	PRIVATE
	E:/VisualStudioFiles/vcpkg/ku/armadillo-14.4.1/lib_win64/libopenblas.lib
)


  # this is heuristically generated, and may not be correct
  find_package(FFTW3 CONFIG REQUIRED)
  target_link_libraries(${PROJECT_NAME} PRIVATE FFTW3::fftw3)

  find_package(FFTW3f CONFIG REQUIRED)
  target_link_libraries(${PROJECT_NAME} PRIVATE FFTW3::fftw3f)

  find_package(FFTW3l CONFIG REQUIRED)
  target_link_libraries(${PROJECT_NAME} PRIVATE FFTW3::fftw3l)
cmake_minimum_required(VERSION 3.14)
set(ROOT_PATH .)
project(WellTrajectoryOptimization CXX)

add_subdirectory(extern/eigen)
add_subdirectory(extern/json)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(./include)
include_directories(./src)

set(HEADERS
	./include/CostFuncs.h 
	./include/CurveHold.h 
	./include/CurveHoldCurveHold.h 
	./include/Hold.h 
	./include/PSO.h
	./include/TrajectoryTemplate.h
	./include/Curve.h
	./include/OptimizeWells.h
	./include/Penalties.h
	./include/Slotting.h
	./include/TestWells.h
	./include/trajectory_test.h
	./include/optimize_wells_test.h
)
	
set(SOURCES
	./src/CostFuncs.cpp 
	./src/CurveHold.cpp 
	./src/CurveHoldCurveHold.cpp  
	./src/Hold.cpp 
	./src/PSO.cpp 
	./src/TrajectoryTemplate.cpp 
	./src/Curve.cpp
	./src/OptimizeWells.cpp
	./src/Penalties.cpp
	./src/Slotting.cpp
	./src/Source.cpp
	./src/TestWells.cpp
	./src/trajectory_test.cpp
	./src/optimize_wells_test.cpp
)
include_directories(${Eigen_BINARY_DIR})
include_directories(${NLOHMANN_JSON_INCLUDE_BUILD_DIR})
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_link_libraries(${PROJECT_NAME} PRIVATE 
		eigen
		nlohmann_json)
cmake_minimum_required(VERSION 3.14)
set(ROOT_PATH .)
project(WellTrajectoryOptimization CXX)
include(FetchContent)
FetchContent_Declare(eigen
	GIT_REPOSITORY "https://gitlab.com/libeigen/eigen"
	GIT_TAG "3bb6a48d8c171cf20b5f8e48bfb4e424fbd4f79e")

FetchContent_MakeAvailable(eigen)

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
	./include/WellTrajectorySolver.h
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
	./src/WellTrajectorySolver.cpp 
	./src/Curve.cpp
	./src/OptimizeWells.cpp
	./src/Penalties.cpp
	./src/Slotting.cpp
	./src/Source.cpp
	./src/TestWells.cpp
	./src/trajectory_test.cpp
	./src/optimize_wells_test.cpp
)
include_directories(${eigen_SOURCE_DIR})
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_link_libraries(${PROJECT_NAME} PUBLIC eigen)

CXX=h5c++

# Optimization flags 
# No optimization, usefull for debugging.
CXX_FLAGS=-std=c++11 -g3 -O0 -Wall
# Optimize for current processor.
# CXX_FLAGS=-std=c++11 -O3

all: 3D.cpp
	$(CXX) $< -o 3D $(CXX_FLAGS)

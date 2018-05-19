CXX=h5c++

# Optimization flags 
# No optimization, usefull for debugging.
CXX_FLAGS=-std=c++11 -g3 -O0 -Wall
# Optimize for current processor.
# CXX_FLAGS=-std=c++11 -O3

all:
	g++ FDTD.cpp -o FDTD -I/usr/include/hdf5/serial -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -Wdate-time -D_FORTIFY_SOURCE=2 -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -std=c++11

luis: 3D.cpp
	$(CXX) $< -o 3D $(CXX_FLAGS)

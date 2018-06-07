CXX=h5c++

# Optimization flags 
# No optimization, usefull for debugging.
CXX_FLAGS=-std=c++11 -g3 -O0 -Wall
# Optimize for current processor.
# CXX_FLAGS=-std=c++11 -O3

all:
	g++ FDTD.cpp -o FDTD -I/usr/include/hdf5/serial -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -Wdate-time -D_FORTIFY_SOURCE=2 -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -std=c++11

luis: FDTD.cpp
	$(CXX) $< -o 3D $(CXX_FLAGS)

plot:
	python plotcm.py Ex_teor XY
	python plotcm.py Ex_teor XZ
	python plotcm.py Ex_teor YZ
	python plotcm.py Ex_sim XY
	python plotcm.py Ex_sim XZ
	python plotcm.py Ex_sim YZ
	python plotcm.py Ey_teor XY
	python plotcm.py Ey_teor XZ
	python plotcm.py Ey_teor YZ
	python plotcm.py Ey_sim XY
	python plotcm.py Ey_sim XZ
	python plotcm.py Ey_sim YZ
	python plotcm.py Ez_teor XY
	python plotcm.py Ez_teor XZ
	python plotcm.py Ez_teor YZ
	python plotcm.py Ez_sim XY
	python plotcm.py Ez_sim XZ
	python plotcm.py Ez_sim YZ

ploterror:
	#python ploterror.py P XY
	#python ploterror.py P XZ
	#python ploterror.py P YZ
	#python ploterror.py Ex XY
	#python ploterror.py Ex XZ
	#python ploterror.py Ex YZ
	#python ploterror.py Ey XY
	#python ploterror.py Ey XZ
	#python ploterror.py Ey YZ
	python ploterror.py Ez XY
	python ploterror.py Ez XZ
	python ploterror.py Ez YZ

fieldlines:
	python field_lines.py E sim 1
	python field_lines.py E sim 1 above
	python field_lines.py E sim 60
	python field_lines.py E sim 60 above
	python field_lines.py E sim 120
	python field_lines.py E sim 120 above
	python field_lines.py H sim 1
	python field_lines.py H sim 1 above
	python field_lines.py H sim 60
	python field_lines.py H sim 60 above
	python field_lines.py H sim 120
	python field_lines.py H sim 120 above
	python field_lines.py P sim 1
	python field_lines.py P sim 1 above
	python field_lines.py P sim 60
	python field_lines.py P sim 60 above
	python field_lines.py P sim 120
	python field_lines.py P sim 120 above
	python field_lines.py EHP sim 1
	python field_lines.py EHP sim 1 above
	python field_lines.py EHP sim 60
	python field_lines.py EHP sim 60 above
	python field_lines.py EHP sim 120
	python field_lines.py EHP sim 120 above
	python field_lines.py E teor 1
	python field_lines.py E teor 1 above
	python field_lines.py E teor 60
	python field_lines.py E teor 60 above
	python field_lines.py E teor 120
	python field_lines.py E teor 120 above
	python field_lines.py H teor 1
	python field_lines.py H teor 1 above
	python field_lines.py H teor 60
	python field_lines.py H teor 60 above
	python field_lines.py H teor 120
	python field_lines.py H teor 120 above
	python field_lines.py P teor 1
	python field_lines.py P teor 1 above
	python field_lines.py P teor 60
	python field_lines.py P teor 60 above
	python field_lines.py P teor 120
	python field_lines.py P teor 120 above
	python field_lines.py EHP teor 1
	python field_lines.py EHP teor 1 above
	python field_lines.py EHP teor 60
	python field_lines.py EHP teor 60 above
	python field_lines.py EHP teor 120
	python field_lines.py EHP teor 120 above

above:
	python above.py E sim 1
	python above.py E sim 60
	python above.py E sim 120
	python above.py H sim 1
	python above.py H sim 60
	python above.py H sim 120
	python above.py P sim 1
	python above.py P sim 60
	python above.py P sim 120
	python above.py EHP sim 1
	python above.py EHP sim 60
	python above.py EHP sim 120
	python above.py E teor 1
	python above.py E teor 60
	python above.py E teor 120
	python above.py H teor 1
	python above.py H teor 60
	python above.py H teor 120
	python above.py P teor 1
	python above.py P teor 60
	python above.py P teor 120
	python above.py EHP teor 1
	python above.py EHP teor 60
	python above.py EHP teor 120
shut:
	shutdown now
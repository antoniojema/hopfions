set xrange[0:100]
set cbrange[0:1]
set key off

set pm3d map
do for [i=0:100]{
	set yrange[i*101:(i+1)*101-1]
	splot "3Dcm.txt" matrix
}

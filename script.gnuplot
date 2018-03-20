set xrange[0:100]
set yrange[-3:3]
set key off

do for [i=0:100]{
	plot "3D.txt" u 1:2 every ::i*101::(i+1)*101-1 w lines
	pause 0.1
}
	

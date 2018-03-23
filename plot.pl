do for [i=0:100] {
    plot "3Dcpp.txt" using 1:2 every ::i*101::(i+1)*101-1 with lines
    pause 0.2
}


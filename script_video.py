import os
import sys
import numpy as np
ARG = sys.argv[1]

for i in np.arange(1,201,10):
	print 'python field_lines_'+ARG+'.py P sim '+str(i)
	os.system('python field_lines_'+ARG+'.py P sim '+str(i))
	print 'python field_lines_'+ARG+'.py P teor '+str(i)
	os.system('python field_lines_'+ARG+'.py P teor '+str(i))
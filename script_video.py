import os
import sys
import numpy as np
import time
ARG = sys.argv[1]
try:
	ARG2 = sys.argv[2]
except:
	ARG2 = None

INTERVAL = np.arange(1,122,4)

if ARG == 'all':
	for i in INTERVAL:
		time.sleep(0.1)
		print 'python field_lines_11.py P sim '+str(i)+' video'
		os.system('python field_lines_11.py P sim '+str(i)+' video')
		print 'python field_lines_11.py P teor '+str(i)+' video'
		os.system('python field_lines_11.py P teor '+str(i)+' video')
		
		print 'python field_lines_23.py P sim '+str(i)+' video'
		os.system('python field_lines_23.py P sim '+str(i)+' video')
		print 'python field_lines_23.py P teor '+str(i)+' video'
		os.system('python field_lines_23.py P teor '+str(i)+' video')

else:
	for i in INTERVAL:
		time.sleep(0.1)
		print 'python field_lines_'+ARG+'.py P sim '+str(i)+' video'
		os.system('python field_lines_'+ARG+'.py P sim '+str(i)+' video')
		print 'python field_lines_'+ARG+'.py P teor '+str(i)+' video'
		os.system('python field_lines_'+ARG+'.py P teor '+str(i)+' video')

time.sleep(5)

if ARG2 == 'shut':
	os.system('shutdown now')
import numpy as np
from PIL import Image
from scipy import misc
import sys

field = sys.argv[1]
which = sys.argv[2]
n = sys.argv[3]
if n == '1':
    t = '-15'

elif n == '60':
    t = '0'

elif n == '120':
    t = '15'

img = misc.imread('../hopfions_memoria/media/11_flow_t'+t+'_'+field+'_'+which+'_above.png')

a = np.swapaxes(img[81:(len(img)-81)], 0, 1)
img = np.swapaxes(a[81:(len(a)-82)], 0, 1)

misc.imsave('../hopfions_memoria/media/11_flow_t'+t+'_'+field+'_'+which+'_above.png',img)
#misc.imsave('prueba.png',img)

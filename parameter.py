
import numpy as np

aA=float(73.2180683898237)
y_xratio=float(1.68)
bB=y_xratio*aA
cC=float(12.388798)
Xlen=np.ceil(aA*2)
Ylen=np.ceil(bB*2)
Zlen=np.ceil(cC*2)
x_min=0
y_min=0
z_min=0
x_max=Xlen
y_max=Ylen
z_max=Zlen
nel_x=Xlen/2 
nel_y=Ylen/2
nel_z=Zlen/2
nel_x=int(nel_x)
nel_y=int(nel_y)
nel_z=int(nel_z)
x_width=(x_max-x_min)/nel_x 
y_width=(y_max-y_min)/nel_y
z_width=(z_max-z_min)/nel_z


print ('number of element',nel_x, nel_y, nel_z)
print ('len', Xlen,Ylen,Zlen)
print ('max',x_max,y_max,z_max)
print ('width',x_width,y_width,z_width)
print ('rx,ry,rx',aA,bB,cC)

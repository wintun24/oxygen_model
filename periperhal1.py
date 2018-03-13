'''
This code is to calculat the volume of element that are not complete cube inside the ellipsoid. 
elementeSeparation 0=the cube is completly outside the ellipsoid, 1= the cube is
on the border of ellipsoid. 2= the cube is completely inside the ellipsoid. 
If it is 2, the vol is x_width*y_widht*z_width. If it is 1, has to calculated the portion of volume that is inside the 
ellipsoid. vol_corrected stores the real volume data. 
'''



import numpy as np

from parameter import * #x_min,y_min,z_min,x_max,y_max,z_max,x_width,y_width,z_width
print aA, bB, cC
print nel_x,nel_y,nel_z
print x_max,y_max,z_max

#halfx=x_max/2
#halfy=y_max/2
#halfz=z_max/2
#a=x_max/2
#b=y_max/2
#c=z_max/2
a=aA
b=bB
c=cC

#plotting node of mesh
nodeOfelement=np.zeros((8,nel_x*nel_y*nel_z))
E1=0

for K1 in range (1,nel_z+1):
    for J1 in range (1,nel_y+1):
        for I1 in range(1,nel_x+1):
           
            nodeOfelement[0,E1] = I1+(nel_x+1)*(J1-1)+(nel_x+1)*(nel_y+1)*(K1-1);
            nodeOfelement[1,E1] = nodeOfelement[0,E1]+1;
            nodeOfelement[2,E1] = nodeOfelement[0,E1]+nel_x+1;
            nodeOfelement[3,E1] = nodeOfelement[2,E1]+1;
            nodeOfelement[4,E1] = nodeOfelement[0,E1]+(nel_x+1)*(nel_y+1);
            nodeOfelement[5,E1] = nodeOfelement[1,E1]+(nel_x+1)*(nel_y+1);
            nodeOfelement[6,E1] = nodeOfelement[2,E1]+(nel_x+1)*(nel_y+1);
            nodeOfelement[7,E1] = nodeOfelement[3,E1]+(nel_x+1)*(nel_y+1);
            
            E1 = E1+1;

nodeOfelement=nodeOfelement.T

#ploting node for element
#x = np.linspace(x_min,x_max, x_max/x_width+1)
x = np.linspace(x_min,x_max, nel_x+1)
#y = np.linspace(y_min, y_max, y_max/y_width+1)
y = np.linspace(y_min, y_max, nel_y+1)
#z = np.linspace(z_min, z_max, z_max/z_width+1)
z = np.linspace(z_min, z_max, nel_z+1)
nodes = np.vstack(np.meshgrid(y,z,x)).reshape(3,-1)

y=np.array(nodes[0,:])
z=np.array(nodes[1,:])
x=np.array(nodes[2,:])

elementSeparation=np.zeros((nel_x*nel_y*nel_z,1))
vol_corrected=np.zeros((nel_x*nel_y*nel_z,1))

for ii in range(0,len(nodeOfelement)):
    
       if ((x[nodeOfelement[ii,0]-1]-a)**2/a**2+(y[nodeOfelement[ii,0]-1]-b)**2/b**2+(z[nodeOfelement[ii,0]-1]-c)**2/c**2)<=1 and ((x[nodeOfelement[ii,1]-1]-a)**2/a**2+(y[nodeOfelement[ii,1]-1]-b)**2/b**2+(z[nodeOfelement[ii,1]-1]-c)**2/c**2)<=1 and ((x[nodeOfelement[ii,2]-1]-a)**2/a**2+(y[nodeOfelement[ii,2]-1]-b)**2/b**2+(z[nodeOfelement[ii,2]-1]-c)**2/c**2)<=1 and ((x[nodeOfelement[ii,3]-1]-a)**2/a**2+(y[nodeOfelement[ii,3]-1]-b)**2/b**2+(z[nodeOfelement[ii,3]-1]-c)**2/c**2)<=1 and ((x[nodeOfelement[ii,4]-1]-a)**2/a**2+(y[nodeOfelement[ii,4]-1]-b)**2/b**2+(z[nodeOfelement[ii,4]-1]-c)**2/c**2)<=1 and ((x[nodeOfelement[ii,5]-1]-a)**2/a**2+(y[nodeOfelement[ii,5]-1]-b)**2/b**2+(z[nodeOfelement[ii,5]-1]-c)**2/c**2)<=1 and ((x[nodeOfelement[ii,6]-1]-a)**2/a**2+(y[nodeOfelement[ii,6]-1]-b)**2/b**2+(z[nodeOfelement[ii,6]-1]-c)**2/c**2)<=1 and ((x[nodeOfelement[ii,7]-1]-a)**2/a**2+(y[nodeOfelement[ii,7]-1]-b)**2/b**2+(z[nodeOfelement[ii,7]-1]-c)**2/c**2)<=1 :
         elementSeparation[ii,0]=2
         vol_corrected[ii,0]=x_width*y_width*z_width
         #print ('element is inside')

       elif ((x[nodeOfelement[ii,0]-1]-a)**2/a**2+(y[nodeOfelement[ii,0]-1]-b)**2/b**2+(z[nodeOfelement[ii,0]-1]-c)**2/c**2)>1 and ((x[nodeOfelement[ii,1]-1]-a)**2/a**2+(y[nodeOfelement[ii,1]-1]-b)**2/b**2+(z[nodeOfelement[ii,1]-1]-c)**2/c**2)>1 and ((x[nodeOfelement[ii,2]-1]-a)**2/a**2+(y[nodeOfelement[ii,2]-1]-b)**2/b**2+(z[nodeOfelement[ii,2]-1]-c)**2/c**2)>1 and ((x[nodeOfelement[ii,3]-1]-a)**2/a**2+(y[nodeOfelement[ii,3]-1]-b)**2/b**2+(z[nodeOfelement[ii,3]-1]-c)**2/c**2)>1 and ((x[nodeOfelement[ii,4]-1]-a)**2/a**2+(y[nodeOfelement[ii,4]-1]-b)**2/b**2+(z[nodeOfelement[ii,4]-1]-c)**2/c**2)>1 and ((x[nodeOfelement[ii,5]-1]-a)**2/a**2+(y[nodeOfelement[ii,5]-1]-b)**2/b**2+(z[nodeOfelement[ii,5]-1]-c)**2/c**2)>1 and ((x[nodeOfelement[ii,6]-1]-a)**2/a**2+(y[nodeOfelement[ii,6]-1]-b)**2/b**2+(z[nodeOfelement[ii,6]-1]-c)**2/c**2)>1 and ((x[nodeOfelement[ii,7]-1]-a)**2/a**2+(y[nodeOfelement[ii,7]-1]-b)**2/b**2+(z[nodeOfelement[ii,7]-1]-c)**2/c**2)>1 :

         #print ('element is outside_no branch_not placenta')
         vol_corrected[ii,0]=0
       else: 

         #print ('element on border')
                 
         elementSeparation[ii,0]=1

         startx=x[nodeOfelement[ii,0]-1]
         endx=x[nodeOfelement[ii,1]-1]

         starty=y[nodeOfelement[ii,0]-1]
         endy=y[nodeOfelement[ii,2]-1]

         startz=z[nodeOfelement[ii,0]-1]
         endz=z[nodeOfelement[ii,4]-1]

         xVector = np.linspace(startx, endx, 5)
         yVector = np.linspace(starty, endy, 5)
         zVector = np.linspace(startz, endz, 5)
         
         nodes = np.vstack(np.meshgrid(xVector,yVector,zVector)).reshape(3,-1).T
         pointcount=0
         for jj in range (0, len(nodes)):
             if (nodes[jj,0]-a)**2/a**2+(nodes[jj,1]-b)**2/b**2+(nodes[jj,2]-c)**2/c**2<=1:
                pointcount=pointcount+1
                                        
 
         vol_corrected[ii,0]=float(pointcount)/float(len(nodes))*(x_width*y_width*z_width)
         


    
'''
f=open("BorderElement.exnode",'w')# this is k value in each node of mesh

f.write (" Group name: Test\n")
f.write (" #Fields=1\n")
f.write (" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
f.write (" x.  Value index=1, #Derivatives=0\n")
f.write (" y.  Value index=2, #Derivatives=0\n")
f.write (" z.  Value index=3, #Derivatives=0\n")
bb=0

for ee in range(0,len(elementSeparation)):
  for ff in range(0,8):
    if elementSeparation[ee,0]==1:
     f.write("Node:  "        "%s\n" %(bb+1))
     f.write("          %s\n" %x[nodeOfelement[ee,ff]-1])
     f.write("          %s\n" %y[nodeOfelement[ee,ff]-1])
     f.write("          %s\n" %z[nodeOfelement[ee,ff]-1])
     bb=bb+1




'''


for c in range(0,len(vol_corrected)):
     if elementSeparation[c,0]==0 and vol_corrected[c,0]!=0:
       print 'ERROR1'
     elif elementSeparation[c,0]==1 and vol_corrected[c,0]==0:
       print 'ERROR2'

     elif elementSeparation[c,0]==2 and vol_corrected[c,0]!=x_width*y_width*z_width:
       print 'ERROR3'
   

f=open("Corrected_Vol.txt",'w')
for ee in range(0,len(vol_corrected)):
  
    f.write("%s %s \n" %(elementSeparation[ee,0],vol_corrected[ee,0]))


f.close()



print np.sum(vol_corrected[:,0])










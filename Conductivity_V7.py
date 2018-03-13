import numpy as np
import math
import numpy.matlib 
import sys
from parameter import *
import time



tic = time.time()

########################DATA READ IN######################################


Branch_El = np.loadtxt('NodeOfBrelement_whole.txt')# List of branch element
Branch_Node=np.loadtxt('CoorOfBrelement_whole.txt')# List of Coor of nodes 
Radius_Br=np.loadtxt('Radius_whole.txt')#List of radius of each element
term_block=np.loadtxt('terminal_block_whole.txt')#block element number where terminal br are located
Cor_Vol=np.loadtxt('Corrected_Vol.txt')# corrected vol at the edge of ellipsoid
#master_vol_voxel=np.loadtxt('master_vol_voxel.txt')
pi=math.pi


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


f=open("Node.txt",'w')


for ee in range(0,len(nodeOfelement)):
  
    f.write("%s %s %s %s %s %s %s %s\n" %(nodeOfelement[ee,0],nodeOfelement[ee,1],nodeOfelement[ee,2],nodeOfelement[ee,3],nodeOfelement[ee,4],nodeOfelement[ee,5],nodeOfelement[ee,6],nodeOfelement[ee,7]))


f.close()

#########################################################################################
wt_diam = np.zeros((nel_x*nel_y*nel_z,1))
br_counting=np.zeros((nel_x*nel_y*nel_z,1))
branchvolume=np.zeros((len(Branch_El),1))
#BLOCK ABOVE
master_vol_voxel=np.zeros((nel_x*nel_y*nel_z, 2))
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Br_count=0
for SS in range (0,len(Branch_El)):

  NODEONE=Branch_Node[Branch_El[SS,0]-1]
  NODETWO=Branch_Node[Branch_El[SS,1]-1]
  
  
  Radius=Radius_Br[SS]
  
  X1=NODEONE
  X2=NODETWO

  r=Radius
  interval=0.1

  points=5
  unit_Vx=[1, 0, 0]


  angle_X1X2 = np.arccos(np.dot(unit_Vx,np.subtract(X2,X1))/(np.linalg.norm(unit_Vx)*np.linalg.norm(np.subtract(X2,X1))) )*180/pi


  axis_rot = np.cross([1,0, 0],np.subtract(X2,X1) )


  length_cyl=np.linalg.norm(np.subtract(X2,X1))
  branch_vol = pi*r**2*length_cyl
  
  branchvolume[SS,0]=branch_vol
 


  x1=np.linspace(interval,np.subtract(length_cyl,interval),10)
  
  if x1[-1]== x1[-2]:
     x1_temp = x1[0:len(x1)-1]
     x1=[]
     x1=x1_temp
     print("CHECK x1")

  if len(x1)==0:
   print("x1 is empty")
 

  yp=np.linspace(-r,r,points)
  zp=np.linspace(-r,r,points)
  [yd,zd]=np.meshgrid(yp,zp)

  counter=0
  xyz1=np.zeros((2,len(yd)*len(yd[0])))
 
  for a in range(0,len(yd[0])):
     for b in range(0,len(yd[0])):
       if r-(math.sqrt(yd[a,b]**2+zd[a,b]**2)) >= 0: 
           
           xyz1[0,counter]=yd[a,b];
           xyz1[1,counter]=zd[a,b];
           counter=counter+1

  xyz1=xyz1[:,0:counter]
 
  xyz=np.matlib.repmat(x1[0], 1, counter)
 
  for c in range (1,len(x1)):
     xyz=np.hstack([xyz, np.matlib.repmat(x1[c],1,counter)])

  xyz1=np.matlib.repmat(xyz1[:],1,len(x1))
 
  xyz=np.vstack((xyz,xyz1))
 
 
  if angle_X1X2!=0 or angle_X1X2!=180:
    axis_rot=axis_rot[:]/np.linalg.norm(axis_rot)
    
    R=np.eye(3)
    axis_rot=axis_rot[:]/np.linalg.norm(axis_rot)
    x=angle_X1X2
   
    for ii in range(0,3):
      v=R[:,ii]
      R[:,ii] = v*math.cos(math.radians(x)) + np.cross(axis_rot,v)*math.sin(math.radians(x)) + sum((axis_rot*v))*(1-math.cos(math.radians(x)))*axis_rot
 
  Rxyzold=np.dot(R,xyz)
  X1=np.array([X1]).T
  xyznew=Rxyzold+X1
  X2=np.array([X2])
 
  if angle_X1X2==0: 
    print('ZERO--')
    xyznew = xyz
   
    xyznew[0,:] = -1*xyznew[0,:]+X2[0,0]
    xyznew[1,:] = xyznew[1,:]+X2[0,1]
    xyznew[2,:] = xyznew[2,:]+X2[0,2]
    
  if angle_X1X2==180:
     print('ONE EIGHTY')
     xyznew = xyz
     xyznew[0,:] = xyz[0,:]+X2[0,0]
     xyznew[1,:]= xyznew[1,:]+X2[0,1]
     xyznew[2,:] = xyznew[2,:]+X2[0,2]
    


  xyznew = xyznew.T

  datapoints = xyznew
  '''
  f=open("datapoints.txt",'w')# this is k in each mesh element

  for x in range(0,len(datapoints)):
    f.write("%0.6f %0.6f %0.6f\n" %(datapoints[x,0],datapoints[x,1],datapoints[x,2]))
  f.close() 
  '''

  Dx=[]
  Dy=[]
  Dz=[]
  for cc in range (0,len(datapoints)):

    if (datapoints[cc,0]-aA)**2/aA**2+(datapoints[cc,1]-bB)**2/bB**2+(datapoints[cc,2]-cC)**2/cC**2<=1:

    

      Dx.append(datapoints[cc][0])
      Dy.append(datapoints[cc][1])
      Dz.append(datapoints[cc][2])
      


  Dx=np.array((Dx))
  Dy=np.array((Dy))
  Dz=np.array((Dz))
  DP=np.vstack((Dx,Dy,Dz)).T

  if len(DP)==0:
     sys.exit("CHECK ERROR")
  if len(DP)>len(datapoints):
   print ('DP error')
  data_in_voxel = np.zeros((nel_x*nel_y*nel_z,1))
  
  vol_voxel = np.zeros((nel_x*nel_y*nel_z,1))
  tem_br_counting=np.zeros((nel_x*nel_y*nel_z,1))
  for i in range(0,len(DP)):
  


   num_x = np.floor((DP[i,0]-x_min)/x_width)+1;
   num_y = np.floor((DP[i,1]-y_min)/y_width)+1;
   num_z = np.floor((DP[i,2]-z_min)/z_width)+1;
   
   if num_x > nel_x:
      num_x = nel_x

   if num_y >= nel_y:
      num_y = nel_y
   
   if num_z >= nel_z:
      num_z = nel_z

   e = ((num_z-1)*nel_x*nel_y + (num_y-1)*nel_x + num_x)-1
   
   data_in_voxel[e,0] = data_in_voxel[e,0]+1 
   
  tot = 0;
  for i in range (0,nel_x*nel_y*nel_z):
    tot = tot + data_in_voxel[i,0]
    

  if tot != len(DP):
    print('Check voxel_distribution - not all datapoints are allocated to voxels')
  

  vol_voxel = data_in_voxel/len(DP)*branch_vol
  br_search=np.where(vol_voxel[:,0]!=0)
  tem_br_counting[br_search]=1
  br_counting[:,0]=br_counting[:,0] + tem_br_counting[:,0]


  master_vol_voxel[:,0] = master_vol_voxel[:,0] + vol_voxel[:,0]
 
  master_vol_voxel[:,1] = master_vol_voxel[:,1] + vol_voxel[:,0]*r*2;
  Br_count=Br_count+1
  print 'Branch Count Done', Br_count

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


f=open("brvol.txt",'w')

for x in range(0,len(branchvolume)):
    f.write("%0.6f \n" %branchvolume[x,0])
f.close() 

f=open("branchperelement.txt",'w')

for x in range(0,len(br_counting)):
    f.write("%0.6f \n" %br_counting[x,0])
f.close() 

f=open("master_vol_voxel.txt",'w')

for x in range(0,len(master_vol_voxel)):
    f.write("%0.6f %0.6f\n" %(master_vol_voxel[x,0],master_vol_voxel[x,1]))
f.close() 
#BLOCK BELOW

x = np.linspace(x_min,x_max, nel_x+1)
y = np.linspace(y_min, y_max, nel_y+1)
z = np.linspace(z_min, z_max, nel_z+1)
nodes = np.vstack(np.meshgrid(y,z,x)).reshape(3,-1)

y=np.array(nodes[0,:])
z=np.array(nodes[1,:])
x=np.array(nodes[2,:])


f=open("xyz.txt",'w')


for ee in range(0,len(x)):
    f.write(" %s   %s       %s\n" %(x[ee],y[ee],z[ee]))
    
f.close()





for i in range (0,len(master_vol_voxel)):
    if master_vol_voxel[i,0]>(x_width*y_width*z_width) and Cor_Vol[i,0]==2:
         master_vol_voxel[i,0] = x_width*y_width*z_width
    elif master_vol_voxel[i,0]>Cor_Vol[i,1] and Cor_Vol[i,0]==1 :
         master_vol_voxel[i,0]=Cor_Vol[i,1]




vol_ratio=np.zeros((len(master_vol_voxel),1))

for VV in range (0,len(master_vol_voxel)):
    vol_ratio[VV,0] = master_vol_voxel[VV,0]/Cor_Vol[VV,1] if Cor_Vol[VV,1]!=0 else 0
    
    if Cor_Vol[VV,1]==0 and master_vol_voxel[VV,0]!=0:
       print ('Cor Vol ERROR',VV)
       #sys.exit(VV)
for i in range (0,nel_x*nel_y*nel_z):
    wt_diam[i,0] = master_vol_voxel[i,1]/master_vol_voxel[i,0] if master_vol_voxel[i,0]!=0 else 0

    if master_vol_voxel[i,0]==0 and master_vol_voxel[i,1]!=0:
       print ('wt_diam ERROR',i)



for Q in range (0,nel_x*nel_y*nel_z):
   if term_block[Q] == 1:
      vol_ratio[Q,0] = 0.4#############P################
      wt_diam[Q,0] = 0.05#########################P###################
'''
f=open("vr_before.txt",'w')


for ee in range(0,len(vol_ratio)):
    f.write(" %s %s  \n" %(vol_ratio[ee,0],wt_diam[ee,0]))
    
f.close()
'''
nodeVR = np.zeros(((nel_x+1)*(nel_y+1)*(nel_z+1), 2))

for i in range(0,len(nodeOfelement)):
   for j in range(0, len(nodeOfelement[i])):
       nodeVR[nodeOfelement[i,j]-1,1]=nodeVR[nodeOfelement[i,j]-1,1]+1;
       nodeVR[nodeOfelement[i,j]-1,0]=(nodeVR[nodeOfelement[i,j]-1,0]*(nodeVR[nodeOfelement[i,j]-1,1]-1)+vol_ratio[i,0])/nodeVR[nodeOfelement[i,j]-1,1]

nodeWD = np.zeros(((nel_x+1)*(nel_y+1)*(nel_z+1), 2))

for i in range(0,len(nodeOfelement)):
   for j in range(0, len(nodeOfelement[i])):
       nodeWD[nodeOfelement[i,j]-1,1]=nodeWD[nodeOfelement[i,j]-1,1]+1;
       nodeWD[nodeOfelement[i,j]-1,0]=(nodeWD[nodeOfelement[i,j]-1,0]*(nodeWD[nodeOfelement[i,j]-1,1]-1)+wt_diam[i,0])/nodeWD[nodeOfelement[i,j]-1,1]



eta1 = 0.5
eta2 = 0.5
eta3 = 0.5

refined_IEN=nodeOfelement.T

BASIC = np.array( [(1-eta1)*(1-eta2)*(1-eta3),   eta1*(1-eta2)*(1-eta3),  (1-eta1)*eta2*(1-eta3),  eta1*eta2*(1-eta3),  (1-eta1)*(1-eta2)*eta3 , eta1*(1-eta2)*eta3	,(1-eta1)*eta2*eta3 , eta1*eta2*eta3])

VR_meshele=np.zeros((nel_x*nel_y*nel_z,1))

for E in range(0,nel_x*nel_y*nel_z): 
    
    VR_FATE=np.array([[nodeVR[refined_IEN[0,E]-1,0]],[nodeVR[refined_IEN[1,E]-1,0] ], [nodeVR[refined_IEN[2,E]-1,0]], [nodeVR[refined_IEN[3,E]-1,0]],[nodeVR[refined_IEN[4,E]-1,0]], [nodeVR[refined_IEN[5,E]-1,0]] , [nodeVR[refined_IEN[6,E]-1,0]], [nodeVR[refined_IEN[7,E]-1,0]]])
        
    VR_FATE=VR_FATE.T
    VR_meshele[E,0] = np.sum(BASIC*VR_FATE) 



WD_meshele=np.zeros((nel_x*nel_y*nel_z,1))

for E in range(0,nel_x*nel_y*nel_z): 
    
    WD_FATE=np.array([[nodeWD[refined_IEN[0,E]-1,0]],[nodeWD[refined_IEN[1,E]-1,0] ], [nodeWD[refined_IEN[2,E]-1,0]], [nodeWD[refined_IEN[3,E]-1,0]],[nodeWD[refined_IEN[4,E]-1,0]], [nodeWD[refined_IEN[5,E]-1,0]] , [nodeWD[refined_IEN[6,E]-1,0]], [nodeWD[refined_IEN[7,E]-1,0]]])
        
    WD_FATE=WD_FATE.T
    WD_meshele[E,0] = np.sum(BASIC*WD_FATE) 
'''
f=open("vr_af.txt",'w')


for ee in range(0,len(VR_meshele)):
    f.write(" %s  %s \n" %(VR_meshele[ee,0],WD_meshele[ee,0]))
    
f.close()
'''
kFINAL=((WD_meshele**2)*(1-VR_meshele)**3)/(180*(VR_meshele)**2)


for i in range(0,len(kFINAL)):
   if kFINAL[i]>0.0519*10 or math.isnan(kFINAL[i]):
      kFINAL[i]=0.0519*10



f=open("ConductivityEachEl_1stsmooth.txt",'w')



for LL in range(0,len(kFINAL)):
    
    f.write("%s\n" %kFINAL[LL,0])

f.close()

nodeK = np.zeros(((nel_x+1)*(nel_y+1)*(nel_z+1), 2))




for i in range(0,len(nodeOfelement)):
   for j in range(0, len(nodeOfelement[i])):
       nodeK[nodeOfelement[i,j]-1,1]=nodeK[nodeOfelement[i,j]-1,1]+1;
       nodeK[nodeOfelement[i,j]-1,0]=(nodeK[nodeOfelement[i,j]-1,0]*(nodeK[nodeOfelement[i,j]-1,1]-1)+kFINAL[i])/nodeK[nodeOfelement[i,j]-1,1]
      



f=open("ConductivityEachNode.exnode",'w')

f.write (" Group name: Conductivity\n")
f.write (" #Fields=2\n")
f.write (" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
f.write (" x.  Value index=1, #Derivatives=0\n")
f.write (" y.  Value index=2, #Derivatives=0\n")
f.write (" z.  Value index=3, #Derivatives=0\n")
f.write (" 2) general, field, rectangular cartesian, #Components=1\n")
f.write ("  1.  Value index=4, #Derivatives=0\n")

for ee in range(0,len(nodeK)):
    f.write("Node:  "        "%s\n" %(ee+1))
    f.write("          %s\n" %x[ee])
    f.write("          %s\n" %y[ee])
    f.write("          %s\n" %z[ee])
    f.write("          %s\n" %nodeK[ee,0])




f.close()

######################################################################################################################################################################

eta1 = 0.5
eta2 = 0.5
eta3 = 0.5

refined_IEN=nodeOfelement.T

BASIC = np.array( [(1-eta1)*(1-eta2)*(1-eta3),   eta1*(1-eta2)*(1-eta3),  (1-eta1)*eta2*(1-eta3),  eta1*eta2*(1-eta3),  (1-eta1)*(1-eta2)*eta3 , eta1*(1-eta2)*eta3	,(1-eta1)*eta2*eta3 , eta1*eta2*eta3])

k_meshele=np.zeros((nel_x*nel_y*nel_z,1))

for E in range(0,nel_x*nel_y*nel_z): 
    
    K_FATE=np.array([[nodeK[refined_IEN[0,E]-1,0]],[nodeK[refined_IEN[1,E]-1,0] ], [nodeK[refined_IEN[2,E]-1,0]], [nodeK[refined_IEN[3,E]-1,0]],[nodeK[refined_IEN[4,E]-1,0]], [nodeK[refined_IEN[5,E]-1,0]] , [nodeK[refined_IEN[6,E]-1,0]], [nodeK[refined_IEN[7,E]-1,0]]])
        
    K_FATE=K_FATE.T
    k_meshele[E,0] = np.sum(BASIC*K_FATE) 



f=open("ConductivityEachEl_2ndsmooth.txt",'w')

for LL in range(0,len(k_meshele)):
    
    f.write("%s\n" %k_meshele[LL,0])

f.close()

######################porosity
porosity=np.zeros((len(master_vol_voxel),1))

for pp in range(0,len(porosity)):
 porosity[pp,0]= 1- VR_meshele[pp,0]
 if math.isnan(VR_meshele[pp,0]):
  print ('VR_meshele ERROR')

f=open("porosityEachEl.txt",'w')

for LL in range(0,len(porosity)):
    
    f.write("%s\n" %porosity[LL,0])

f.close()

nodeP = np.zeros(((nel_x+1)*(nel_y+1)*(nel_z+1), 2))

for i in range(0,len(nodeOfelement)):
   for j in range(0, len(nodeOfelement[i])):
       nodeP[nodeOfelement[i,j]-1,1]=nodeP[nodeOfelement[i,j]-1,1]+1;
       nodeP[nodeOfelement[i,j]-1,0]=(nodeP[nodeOfelement[i,j]-1,0]*(nodeP[nodeOfelement[i,j]-1,1]-1)+porosity[i,0])/nodeP[nodeOfelement[i,j]-1,1]

f=open("PorosityEachNode.exnode",'w')

f.write (" Group name: Porosity\n")
f.write (" #Fields=2\n")
f.write (" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
f.write (" x.  Value index=1, #Derivatives=0\n")
f.write (" y.  Value index=2, #Derivatives=0\n")
f.write (" z.  Value index=3, #Derivatives=0\n")
f.write (" 2) general, field, rectangular cartesian, #Components=1\n")
f.write ("  1.  Value index=4, #Derivatives=0\n")

for ee in range(0,len(nodeP)):
    f.write("Node:  "        "%s\n" %(ee+1))
    f.write("          %s\n" %x[ee])
    f.write("          %s\n" %y[ee])
    f.write("          %s\n" %z[ee])
    f.write("          %s\n" %nodeP[ee,0])




f.close()


p_meshele=np.zeros((nel_x*nel_y*nel_z,1))

for E in range(0,nel_x*nel_y*nel_z): 
    
    P_FATE=np.array([[nodeP[refined_IEN[0,E]-1,0]],[nodeP[refined_IEN[1,E]-1,0] ], [nodeP[refined_IEN[2,E]-1,0]], [nodeP[refined_IEN[3,E]-1,0]],[nodeP[refined_IEN[4,E]-1,0]], [nodeP[refined_IEN[5,E]-1,0]] , [nodeP[refined_IEN[6,E]-1,0]], [nodeP[refined_IEN[7,E]-1,0]]])
        
    P_FATE=P_FATE.T
    p_meshele[E,0] = np.sum(BASIC*P_FATE) 



f=open("porosityEachEl_2ndsmooth.txt",'w')

for LL in range(0,len(p_meshele)):
    
    f.write("%s\n" %p_meshele[LL,0])

f.close()



toc = time.time() - tic
print 'elapsed_time',toc


import numpy as np
from parameter import *
'''
x_min=0
y_min=0
z_min=0
x_max=130
y_max=220
z_max=30
x_width=1
y_width=1
z_width=1
nel_x=x_max/x_width
nel_y=y_max/y_width
nel_z=z_max/z_width
'''
Br_El=np.loadtxt('NodeOfBrelement_whole.txt')##########branch element
node=np.loadtxt('CoorOfBrelement_whole.txt')############branch node

#if it is 1, it is terminal br, if it is zero it is not
T_Br=np.ones((len(Br_El),1))
term_block = np.zeros((nel_x*nel_y*nel_z,1))
T_ELIN= np.zeros((nel_x*nel_y*nel_z,1))



for i in range(0, len(Br_El)):
    for j in range(0,len(Br_El)):     
        
        if Br_El[j,0] == Br_El[i,1]:
           T_Br[i] = 0
             
         
Terminal_El_num=np.array(np.where(T_Br[:,0]==1))### this is minus one value


f=open("Terminal_El_num.txt",'w')

for P in range(0,len(Terminal_El_num[0])):
   
    f.write("%d\n"  %(Terminal_El_num[0,P] ))
f.close()
raise SystemExit

#print Terminal_El_num [0,0]


for j in range(0,len(Terminal_El_num[0])):


  Endpoints = node[[Br_El[Terminal_El_num[0][j],1]-1]]
  #print (Endpoints[0,0])
  num_x = np.floor((Endpoints[0,0]-x_min)/x_width)+1;
  num_y = np.floor((Endpoints[0,1]-y_min)/y_width)+1;
  num_z = np.floor((Endpoints[0,2]-z_min)/z_width)+1;
   
  if num_x > nel_x:
     num_x = nel_x

  if num_y >= nel_y:
     num_y = nel_y
   
  if num_z >= nel_z:
     num_z = nel_z

  T_e = ((num_z-1)*nel_x*nel_y + (num_y-1)*nel_x + num_x)-1###############CHECK##########e is numbering of mesh
  term_block[T_e,0] = 1;# if 1, that block element contain terminal branch
  T_ELIN[T_e,0]=Terminal_El_num[0][j]



f=open("terminal_block_whole.txt",'w')

for P in range(0,len(term_block)):
   
    f.write("%d\n"  %(term_block[P] ))
f.close()

         
f=open("terminal_element_whole.txt",'w')# this value is already minus one

for S in range(0,len(Terminal_El_num[0])):
   
    f.write("%d\n"  %(Terminal_El_num[0][S]+1 ))
f.close()
'''

np.nonzero(term_block)
np.nonzero(T_ELIN)
f=open("terminal_element_bl.txt",'w')# this value is already minus one

for O in range(0,len(term_block)):
   
    f.write("%d %d\n"  %(term_block[O], T_ELIN[O]+1 ))
f.close()

'''


  

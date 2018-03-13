
clear all
close all
%% Read in exnode and radius
fid=fopen('vessels.exnode');
data = textscan(fid, '%s', 'delimiter', '\n');
sdata=regexp(data{1,1},' ','split');



%% Parameters
y_xratio=1.68;
counter=0;
last_ven_elem = 127746;%last vein element in vessel info
first_ven_elem  = 63874;

node=zeros(4,1e6);%initialization

for a=1:size(sdata,1)    
    idx1 = find(strcmp(sdata{a}, 'Node:'));
     if isequal(idx1,1)==1  
         counter=counter+1;
         tt=find(~cellfun(@isempty,sdata{a,1}));          
         node(1,counter)=str2double(sdata{a,1}{1,tt(1,2)});%node number
         tt1=find(~cellfun(@isempty,sdata{a+1,1})); %num of version
         
         %disp(tt1)
         
         if(tt1==1) 
           node(2,counter)=str2double(sdata{a+1,1}); %x
           node(3,counter)=str2double(sdata{a+2,1}); %y
           node(4,counter)=str2double(sdata{a+3,1}); %z  
           node(5,counter)=str2double(sdata{a+4,1}{1,1});%radius       
         elseif (tt1(1,1)==1 && tt1(1,2)==4)    %%%%%%%%%%%check      
           node(2,counter)=str2double(sdata{a+1,1}{1,tt(1,1)}); %x
           node(3,counter)=str2double(sdata{a+2,1}{1,tt(1,1)}); %y
           node(4,counter)=str2double(sdata{a+3,1}{1,tt(1,1)});
           node(5,counter)=str2double(sdata{a+4,1}{1,tt(1,1)});%z   
         elseif (tt1(1,1)==1 && tt1(1,2)==3)    %%%%%%%%%%%check  
             
           node(2,counter)=str2double(sdata{a+1,1}{1,tt(1,1)}); %x
           node(3,counter)=str2double(sdata{a+2,1}{1,tt(1,1)}); %y
           node(4,counter)=str2double(sdata{a+3,1}{1,tt(1,1)});
           node(5,counter)=str2double(sdata{a+4,1}{1,tt(1,1)});
         end
     end
end
clearvars data sdata;
nnode = counter;
node(:,node(1,:)==0)=[]; %remove unfilled nodes
node = node';

check=node(:,2:4);

all_zeros = ~any(check, 2);
if sum(all_zeros)>0
    disp("ERROR NODE READ");
    return
end

firsteight=[1;2;3;4;first_ven_elem+1;first_ven_elem+2;first_ven_elem+3;first_ven_elem+4];
check(firsteight,:)=[];
minx=min(check(:,1));
maxx=max(check(:,1));
miny=min(check(:,2));
maxy=max(check(:,2));
minz=min(check(:,3));
maxz=max(check(:,3));
x_length=abs(minx)+maxx;
y_length=abs(miny)+maxy;
z_length=abs(minz)+maxz;
x_radius=x_length/2;
y_radius=y_length/2;
z_radius=z_length/2;

x_extra1=maxx-x_radius;
x_extra2=x_radius-abs(minx);
y_extra1=maxy-y_radius;
y_extra2=y_radius-abs(miny);
z_extra1=maxz-z_radius;
z_extra2=z_radius-abs(minz);

node(:,2)=node(:,2)-x_extra1;
node(:,3)=node(:,3)-y_extra1;
node(:,4)=node(:,4)-z_extra1;
% scatter3(node(:,2),node(:,3),node(:,4),'r.');
%  hold on
% scatter3(check(:,1),check(:,2),check(:,3),'*b');
%view([0 0]);
xr=x_radius;
yr=xr*y_xratio;
zr=max([abs(minz),maxz]);

for p=1:size(check,1)
    xp=node(p,2);
    yp=node(p,3);
    zp=node(p,4);
    
    
   if ( xp^2/xr^2) + ( yp^2/yr^2) + ( zp^2/zr^2) > 1
       
       %disp( ((xp)^2/xr^2)+( (yp)^2/yr^2) +   ( (zp)^2/zr^2))
       outsidepoint(p)=p;
    
   end


end
Ind=find(outsidepoint~=0);

fileID=fopen('nodeeliminated.txt','w');


for aa=1:size(node,1)
    
if ~ismember(aa, Ind(:))
    
fprintf(fileID,'%2.4f %2.4f %2.4f \n',node(aa,2),node(aa,3),node(aa,4));
end

end
fclose(fileID);

  
%%% Read in exelem for connectivity
fid=fopen('vessels.exelem');
e_counter=0; %counter
elem=zeros(3,1e6);%initialise
adata=textscan(fid,'%s','delimiter','\n');
sdata=regexp(adata{1,1},' ','split');
for a=1:size(sdata,1)
    idx1 = find(strcmp(sdata{a}{1,1}, 'Element:'));
    if isequal(idx1,1)==1
       e_counter=e_counter+1;
       tt=find(~cellfun(@isempty,sdata{a,1}));
       elem(1,e_counter)=str2double(sdata{a,1}{1,tt(1,2)}); %elem number
       tt=find(~cellfun(@isempty,sdata{a+2,1}));
       elem(2,e_counter)=str2double(sdata{a+2,1}{1,tt(1,1)}); %node 1
       elem(3,e_counter)=str2double(sdata{a+2,1}{1,tt(1,2)}); %node 2
    end
end
clearvars adata sdata tt;
nel = e_counter;
elem(:,elem(1,:)==0)=[]; %remove unfilled nodes

elem = elem';


fileID=fopen('CoorOfBrelement_whole.txt','w');


for aa=1:size(node,1)
    

    
fprintf(fileID,'%2.4f %2.4f %2.4f \n',node(aa,2)+xr,node(aa,3)+yr,node(aa,4)+zr);%%summation may change


end
fclose(fileID);

fileID=fopen('NodeOfBrelement_whole.txt','w');


for bb=1:last_ven_elem
    
if ~ismember(elem(bb,2), Ind(:)) && ~ismember(elem(bb,3), Ind(:))
    
fprintf(fileID,'%d %d \n',elem(bb,2),elem(bb,3));
elementinclude(bb)=bb;

end
end
fclose(fileID);

 toinclude=(find(elementinclude~=0));


D=textread('radius.exelem','%s','delimiter','\n');

m=1;

for n=11:3:size(D,1)
nodepair=D{n};


X(m,:)= str2num(nodepair);

m=m+1;
end
X=X(:,1);


fileID=fopen('Radius_whole.txt','w');


for cc=1:last_ven_elem
   if ismember(cc, toinclude(:))
   
fprintf(fileID,'%0.4f\n',X(cc,1));
    end

end
fclose(fileID);


e=load('nodeeliminated.txt');
scatter3(e(:,1),e(:,2),e(:,3),'r.');
hold on
[x, y, z] = ellipsoid(0,0,0,xr,yr,zr,50);
surf(x, y, z)
view([0 90]);
pi=3.1415927;
vol_placenta=(4/3)*pi*xr*yr*zr;


Coor=load('CoorOfBrelement_whole.txt');
El=load('NodeOfBrelement_whole.txt');
R_br=load('Radius_whole.txt');
for i= 1:size(El,1)
   Node1=Coor(El(i,1),:);
   N1x=Node1(1,1);
   N1y=Node1(1,2);
   N1z=Node1(1,3);
   Node2=Coor(El(i,2),:);
   N2x=Node2(1,1);
   N2y=Node2(1,2);
   N2z=Node2(1,3);
   
   length=sqrt((N1x-N2x)^2+(N1y-N2y)^2+(N1z-N2z)^2);
   r_br=R_br(i);
   vol(i)=pi*r_br^2*length;
    
    
end
vol=vol';
vascular_vol=sum(vol);


fileID=fopen('V.exnode','w');
 nodenumber=1;
 
 fprintf(fileID,'  Group name: artery\n');
 fprintf(fileID,'  #Fields=2\n');
 fprintf(fileID,'  1) coordinates, coordinate, rectangular cartesian, #Components=3\n');
 fprintf(fileID,'  x.  Value index=1, #Derivatives=0\n');
 fprintf(fileID,'  y.  Value index=2, #Derivatives=0\n');
 fprintf(fileID,'  z.  Value index=3, #Derivatives=0\n');
 fprintf(fileID,' 2) general, field, rectangular cartesian, #Components=1\n');
 fprintf(fileID,'  1.  Value index=4, #Derivatives=0\n');
 
 for i=1:size(node,1)
    % if ~ismember(i, Ind(:))
     
fprintf(fileID,'               Node: %d\n',nodenumber);
     
 fprintf(fileID,'%4.4f \n',node(i,2));
 fprintf(fileID,'%4.4f \n',node(i,3));
 fprintf(fileID,'%4.4f \n',node(i,4));
 fprintf(fileID,'%4.2f \n',node(i,5));
 nodenumber=nodenumber+1;
     end

 
 fclose(fileID);
 
 fileID=fopen('V.exelem','w');
 elnumber=1;
fprintf(fileID,' Group name: artery\n');
fprintf(fileID,' Shape. Dimension=1, line\n');
fprintf(fileID,' #Scale factor sets=1\n');
fprintf(fileID,' l.Lagrange, #Scale factors=2\n');
fprintf(fileID,' #Nodes=2\n');
fprintf(fileID,' #Fields=1\n');
fprintf(fileID,' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n');
fprintf(fileID,' x. l.Lagrange, no modify, standard node based.\n');
fprintf(fileID,'  #Nodes=2\n');
fprintf(fileID,'  1. #Values=1\n');
fprintf(fileID,'    Value indices: 1\n');
fprintf(fileID,'    Scale factor indices: 1\n');
fprintf(fileID,'  2. #Values=1\n');
fprintf(fileID,'    Value indices: 1\n');
fprintf(fileID,'    Scale factor indices: 2\n');
fprintf(fileID,' y. l.Lagrange, no modify, standard node based.\n');
fprintf(fileID,'  #Nodes=2\n');
fprintf(fileID,'  1. #Values=1\n');
fprintf(fileID,'    Value indices: 1\n');
fprintf(fileID,'    Scale factor indices: 1\n');
fprintf(fileID,'  2. #Values=1\n');
fprintf(fileID,'    Value indices: 1\n');
fprintf(fileID,'    Scale factor indices: 2\n');
fprintf(fileID,' z. l.Lagrange, no modify, standard node based.\n');
fprintf(fileID,'  #Nodes=2\n');
fprintf(fileID,'  1. #Values=1\n');
fprintf(fileID,'    Value indices: 1\n');
fprintf(fileID,'    Scale factor indices: 1\n');
fprintf(fileID,'  2. #Values=1\n');
fprintf(fileID,'    Value indices: 1\n');
fprintf(fileID,'    Scale factor indices: 2\n');
 
 
 for bb=1:last_ven_elem
 
 
 
 if ~ismember(elem(bb,2), Ind(:)) && ~ismember(elem(bb,3), Ind(:))
fprintf(fileID,' Element: %d 0 0\n',elnumber);
fprintf(fileID,' Nodes:\n');
fprintf(fileID,' %4d %4d\n',elem(bb,2),elem(bb,3));
fprintf(fileID,' Scale factors:\n');
fprintf(fileID,' 1.000000e+00    1.000000e+00\n');

elnumber=elnumber+1;

 end
 end
 
 fclose(fileID);
 
 
 
 
 
 
 
 
 
 

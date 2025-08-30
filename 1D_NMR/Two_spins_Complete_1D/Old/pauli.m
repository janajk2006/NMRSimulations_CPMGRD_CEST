a=[1 0]';
b=[0 1]'; 

Ii=[1 0 
    0 1]; 
I=Ii; 

Px=[0 1 
    1 0]; 

Py=[ 0 -i 
     i 0]; 

Pz=[1 0 
    0 -1]; 

Pa=[1 0 
    0 0]; 

Pb=[0 0 
    0 1]; 

Ix=Px/2; 
Iy=Py/2; 
Iz=Pz/2; 
Ia=Pa; 
Ib=Pb; 
Ip=Ix+i*Iy; 
Im=Ix-i*Iy; 
H=[1 1 
   1 -1]/sqrt(2);

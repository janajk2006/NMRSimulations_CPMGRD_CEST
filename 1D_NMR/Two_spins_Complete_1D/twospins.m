pauli; 
Iii=diprod(Ii,Ii); 

Iia=diprod(Ii,Ia); 
Iai=diprod(Ia,Ii); 
Ibi=diprod(Ib,Ii); 
Iib=diprod(Ii,Ib);
Iaa=diprod(Ia,Ia); 
Iab=diprod(Ia,Ib); 
Iba=diprod(Ib,Ia); 
Ibb=diprod(Ib,Ib); 
Ixi=diprod(Ix,Ii); 
Iyi=diprod(Iy,Ii); 
Izi=diprod(Iz,Ii); 
Iix=diprod(Ii,Ix); 
Iiy=diprod(Ii,Iy); 
Iiz=diprod(Ii,Iz); 
Ixx=diprod(Ix,Ix); 
Ixy=diprod(Ix,Iy); 
Ixz=diprod(Ix,Iz); 
Iyx=diprod(Iy,Ix); 
Iyy=diprod(Iy,Iy); 
Iyz=diprod(Iy,Iz); 
Izx=diprod(Iz,Ix); 
Izy=diprod(Iz,Iy); 
Izz=diprod(Iz,Iz); 

Ipi=Ixi+i*Iyi; 
Iip=Iix+i*Iiy; 
Imi=Ixi-i*Iyi; 
Iim=Iix-i*Iiy; 

Hhh=diprod(H,H); 
Hih=diprod(Ii,H); 
Hhi=diprod(H,Ii); 


% Single transition operators:
% ---------------------------- 
Ixa=diprod(Ix,Ia); 
Ixb=diprod(Ix,Ib); 
Iax=diprod(Ia,Ix); 
Ibx=diprod(Ib,Ix); 


Iya=diprod(Iy,Ia); 
Iyb=diprod(Iy,Ib); 
Iay=diprod(Ia,Iy); 
Iby=diprod(Ib,Iy); 


Iza=diprod(Iz,Ia); 
Izb=diprod(Iz,Ib); 
Iaz=diprod(Ia,Iz); 
Ibz=diprod(Ib,Iz); 

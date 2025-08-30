
pauli

% Cartesian operators: i,x,y,z

HiNi=diprod(Hi,Ni);
HiNx=diprod(Hi,Nx);
HiNy=diprod(Hi,Ny);
HiNz=diprod(Hi,Nz);
HxNi=diprod(Hx,Ni);
HxNx=diprod(Hx,Nx);
HxNy=diprod(Hx,Ny);
HxNz=diprod(Hx,Nz);
HyNi=diprod(Hy,Ni);
HyNx=diprod(Hy,Nx);
HyNy=diprod(Hy,Ny);
HyNz=diprod(Hy,Nz);
HzNi=diprod(Hz,Ni);
HzNx=diprod(Hz,Nx);
HzNy=diprod(Hz,Ny);
HzNz=diprod(Hz,Nz);

% Ladder operators: i,+,-

HiNi=diprod(Hi,Ni);
HiNp=diprod(Hi,Np);
HiNm=diprod(Hi,Nm);
HpNi=diprod(Hp,Ni);
HpNp=diprod(Hp,Np);
HpNm=diprod(Hp,Nm);
HmNi=diprod(Hm,Ni);
HmNp=diprod(Hm,Np);
HmNm=diprod(Hm,Nm);

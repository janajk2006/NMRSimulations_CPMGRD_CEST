
%Relaxation matrix using Redfield theory

% transcribed by Jana.K
% Ackowledgement:  ScottA. Smith, Tilo Levante, T. S. Mahesh, Ilya Kuprov 
% =======================================================================

function [L_DD_RF]=Relax_matrix(v_ppm,offset_ppm,RF_params,Secular_approx,J_H1N2,n_spins,MHZ_H,v_iso,v_delzz,v_eta,Coord_1,Coord_2,EA_1,EA_2,S,Tau);

nspins


gamB1=RF_params(1,1);
Wrf=RF_params(1,2);



Tau_hsqc=1/(4*J_H1N2);
Dim=power(n_spins,2);


% Calculating the MHZ for H,15N for the desired MHZ with respecto values at 400 MHz

Rfreq=[400.130 -40.561];	% Ref spectrometer 400 MHz H,15N
Gamma_H=2.67515255*power(10,8);
Rfreq_H=400.130;

Rfreq_N=-40.561;
Gamma_N=Rfreq_N*Gamma_H/Rfreq_H;

Gamma_1=Rfreq.*Gamma_H/Rfreq_H;  % its an array, contains info on H, N as first and sencond elements
MHZ=(Gamma_1.*(MHZ_H/Gamma_H)); %  its an array, contains info on H, N as first and sencond elements







% INITIAL DENSITY MATRIX AND HAMILTONIANS 
% --------------------------------------- 

offset_Hz=(diag(offset_ppm*eye(n_spins,n_spins))*MHZ')';
v_Hz=(diag(v_ppm*eye(n_spins,n_spins))*MHZ')';
v_Hz_offseted=v_Hz-offset_Hz;

H_csham=-(v_Hz_offseted(1,1)*HzNi); % Chemical Shift Hamiltonian 
N_csham=-(v_Hz_offseted(1,2)*HiNz); % Chemical Shift Hamiltonian 

T_csham=H_csham+N_csham; % Chemical Shift Hamiltonian 
T_jham=J_H1N2*HzNz;% Scalar coupling Hamiltonian 
T_Jham=J_H1N2*(HxNx+HyNy+HzNz);% Scalar coupling Hamiltonian % use full Jham for CSA relax effect
T_ham=T_csham+T_jham;


% RF related modifications in Hamiltonian

H1rot = gamB1*(HiNx);		% H1 for field along the x-axis
Heff = T_ham - Wrf*(HiNz)+ H1rot;	% Ho to iso rotating frame and Add in the B1 field term







% Liouville parameters
%---------------------

DIM=Dim*Dim;

L_DD_RF=zeros(DIM,DIM);
%L_DC=zeros(DIM,DIM);
%L_CC=zeros(DIM,DIM);

% Relaxation general parameters
%------------------------------

mu = 1.e-7;			% [mu/(4*pi)] J-sec -C  -m
HBAR = 1.05457266e-34;	% Plancks constant (h/2PI)     (J-sec)
mu0d4pi = 1.0e-7;		% mu0/4p (J-sec C  m  )
pi2 = 6.283185307;		% 2*pi
hbar = 1.05459e-34;		% hbar (J-sec)

RAD2HZ=1/(2*pi());
HZ2RAD=2*pi();
DEG2RAD=pi()/180;


v_aniso=v_delzz.*3/2;
Coords= [Coord_1; Coord_2].*power(10,-10);	% In angs
EA= [EA_1;EA_2];


% Calculate Theta, Phi and distance matrices for dipolar mechanism
Thetas=zeros(n_spins,n_spins);
Phis=zeros(n_spins,n_spins);
Distance=zeros(n_spins,n_spins);

X=0;
Y=0;
Z=0;
Rads=0;
for i=1:n_spins
for j=1:n_spins
	Rads=sqrt(power(Coords(i,1)-Coords(j,1),2)+power(Coords(i,2)-Coords(j,2),2)+power(Coords(i,3)-Coords(j,3),2));
	Distance(i,j)=Rads;

	if Rads==0
	Thetas(i,j)=0;
	else
	Thetas(i,j)=acos((Coords(j,3)-Coords(i,3))/Rads);
	end

	%printf("Pt1: x:%0.12lf y:%0.12lf z:%0.12lf\n",Coords(i,1),Coords(i,2),Coords(i,3));
	%printf("Pt2: x:%0.12lf y:%0.12lf z:%0.12lf\n",Coords(j,1),Coords(j,2),Coords(j,3));

	X=Coords(j,1)-Coords(i,1);
	Y=Coords(j,2)-Coords(i,2);

	if Y==0
		if X>=0
		Phis(i,j)=0;
		else
		Phis(i,j)=pi;
		end
	elseif Y<0
		if X==0
		Phis(i,j)=1.5*pi;
		elseif X<0
		Phis(i,j)=atan(Y/X)+pi;	
		else
		Phis(i,j)=atan(Y/X)+2*pi;
		end
	else
		if X==0
		Phis(i,j)=0.5*pi;
		elseif X<0
		Phis(i,j)=atan(Y/X)+pi;	
		else
		Phis(i,j)=atan(Y/X);
		end
	end


end %j
end %i


%-------------------------%
% Dipolar mechanism	  %
%-------------------------%
% as implemented in GAMMA %
%-------------------------%

% Calculate diffusion constants (Used by CSA-DD also): diff_x,y,z;Tau_eff

diff_x=1/(6*Tau(1,1)); % correlation along x,y,z
diff_y=1/(6*Tau(1,2)); % correlation along x,y,z
diff_z=1/(6*Tau(1,3)); % correlation along x,y,z

diff_p=0.5*(diff_x+diff_y);
diff_m=0.5*(diff_x-diff_y);
diff_1=(diff_z-diff_p);
diff_2=2*sqrt(power(diff_1,2)+3.0*power(diff_m,2));
Tau_eff=zeros(1,5);
Tau_eff(1,1)=1/(2*diff_p+4*diff_z);
Tau_eff(1,2)=1/(5*diff_p-3*diff_m+diff_z);
Tau_eff(1,3)=1/(4*diff_p+2*diff_z-diff_2);
Tau_eff(1,4)=1/(5*diff_p+3*diff_m+diff_z);
Tau_eff(1,5)=1/(4*diff_p+2*diff_z+diff_2);



% Caluculate Chi constant (Used by CSA-DD also):
arg1=0;
chi=1.5707963;	%pi/2
if diff_z~=diff_p
arg1=sqrt(3)*diff_m/(diff_z-diff_p);
chi=-atan(arg1);
end



% Calculate dipoles

dip_As=zeros(1,(n_spins-1)*(n_spins-1));	% Dipolar coupling constant of spin pairs
dip_idx=zeros((n_spins-1)*(n_spins-1),2); % index of spin pairs
dip_T=zeros((n_spins-1)*(n_spins-1),5);	% Dipolar spherical tensor spatial 
dip=1;
K = mu*HBAR;

for i=1:n_spins-1
for j=i+1:n_spins

Rads=Distance(i,j);
alpha_temp=Phis(i,j);
beta_temp=Thetas(i,j);
delzz_D=Gamma_1(1,i)*Gamma_1(1,j)*K/power(Rads,3);
iso_D=0;
eta_D=0;

dip_T(dip,1)=0.5*delzz_D*eta_D;	% l=2;m=2
dip_T(dip,2)=0;			% l=2;m=1
dip_T(dip,3)=sqrt(3/2)*delzz_D;	% l=2;m=0
dip_T(dip,4)=0;			% l=2;m=-1
dip_T(dip,5)=0.5*delzz_D*eta_D;	% l=2;m=-2
dip_As(1,dip)=delzz_D;		% Dipolar coupling constant
dip_idx(dip,1)=i;
dip_idx(dip,2)=2;
dip=dip+1;

end
end
dip=dip-1;
dip_T(1:dip,:);
dip_As(1,1:dip);

% Setup HZ

HZ=-MHZ(1,1)*power(10,6)*HzNi-MHZ(1,2)*power(10,6)*HiNz;
%HZ=H_V'*HZ*H_V;
w_omega_D=zeros(1,Dim);

for i=1:Dim
%w_omega_D(1,i)=real(HZ(i,i))+ real(H_D(i,i));
w_omega_D(1,i)=real(HZ(i,i))+ real(T_ham(i,i));
%printf("w_omega_D[%d]:%f (HZ:%f  Ho:%f)\n",i,w_omega_D(1,i),real(HZ(i,i)),real(H_D(i,i)));
end


% Dipolar interction constants between all spin pairs by calculating inter nuclear distances
% required for DD-CSA also

xiD=zeros(n_spins,n_spins);
a=1.941625913;
K = -2.0*a*mu0d4pi*hbar;

for i=1:n_spins
for j=1:n_spins
	if i~=j
	xiD(i,j)=K*Gamma_1(1,i)*Gamma_1(1,j)/(power(Distance(i,j),3));
	end
end
end


% Generate the spin tensor in AAS
% required for DD-CSA also

T12_D_20=zeros(Dim,Dim);
T12_D_2_1=zeros(Dim,Dim);
T12_D_2_2=zeros(Dim,Dim);
T12_D_21=zeros(Dim,Dim);
T12_D_22=zeros(Dim,Dim);


T12_D_20=1/sqrt(6)*(2*HzNi*HiNz-0.5*(HpNi*HiNm+HmNi*HiNp));
T12_D_21=-0.5*(HpNi*HiNz+HzNi*HiNp);
T12_D_2_1=0.5*(HmNi*HiNz+HzNi*HiNm);
T12_D_22=0.5*(HpNi*HiNp);
T12_D_2_2=0.5*(HmNi*HiNm);

%{
T12_D_20
T12_D_21
T12_D_2_1
T12_D_22
T12_D_2_2
pause
%}
T_D=zeros(Dim,Dim,5,n_spins,n_spins);

for i=1:2
for j=1:2
	if i==j
	T_D(:,:,1,i,j)=eye(Dim,Dim);
	T_D(:,:,2,i,j)=eye(Dim,Dim);
	T_D(:,:,3,i,j)=eye(Dim,Dim);
	T_D(:,:,4,i,j)=eye(Dim,Dim);
	T_D(:,:,5,i,j)=eye(Dim,Dim);
	else 
		if i==1 && j==2 
		T_D(:,:,1,i,j)=T12_D_2_2(:,:);
		T_D(:,:,2,i,j)=T12_D_2_1(:,:);
		T_D(:,:,3,i,j)=T12_D_20(:,:);
		T_D(:,:,4,i,j)=T12_D_21(:,:);
		T_D(:,:,5,i,j)=T12_D_22(:,:);
		elseif  i==2 && j==1
		T_D(:,:,1,i,j)=T12_D_2_2(:,:);
		T_D(:,:,2,i,j)=T12_D_2_1(:,:);
		T_D(:,:,3,i,j)=T12_D_20(:,:);
		T_D(:,:,4,i,j)=T12_D_21(:,:);
		T_D(:,:,5,i,j)=T12_D_22(:,:);
		end %if 12,13,23
	end %i==j
end %j
end %i




% Refield matrix element for DD mechanism

xi1=0;
xi2=0;
xi1xi2=0;
alphaij=0;
betaij=0;
alphakl=0;
betakl=0;

ij=1;
kl=0;
w=0;
w0=0;
w1=0;
w2=0;
wi=0;
wj=0;
rank=2;
J=zeros(Dim,Dim);
J_temp=0;
EA1=zeros(1,3);
EA2=zeros(1,3);
c1s=zeros(1,5);
c2s=zeros(1,5);
eta=0;


% For RF field

%[Heff_V,Heff_D]=eig(Heff);
w_omega_rf=zeros(1,Dim);
for i=1:Dim
%w_omega_rf(1,i)=real(Heff_D(i,i));
w_omega_rf(1,i)=real(Heff(i,i));
%printf("w_omega[%d]:%f (HZ:%f  Ho:%f)\n",i,w_omega(1,i),real(HZ(i,i)),real(H_D(i,i)));
end


L_RF=zeros(DIM,DIM);
L_ROT=zeros(DIM,DIM);


mWrf = 0;
J=zeros(Dim,Dim,5);	% for RF effect
J_temp=0;

l3=0;
delwaa=0;
delwbb=0;



for i=1:n_spins-1
for j=i+1:n_spins

	xi1=xiD(i,j);
	alphaij=Phis(i,j);
	betaij=Thetas(i,j);
	gammaij=0;

	% spectral density and spatial spherical harmonics coefficient combined for 5 elements of spatial operator :: five spectral density function coefficients  for isotropic rigid diffusion
	% EA   : Euler angles (radians)
	% chi  : Ratio of rotational correlation times
	% eta  : Tensor value


	sqrt3 = sqrt(3.0);
	K0 = 0.546274215;		% sqrt[15/(16*pi)]
	K1 = 1.092548431;		% sqrt[15/(4.0*pi)]
	K2 = 0.315391565;		% sqrt[5/(16.0*pi)]

	S2g = sin(2.0*gammaij);
	C2g = cos(2.0*gammaij);
	S2a = sin(2.0*alphaij);
	C2a = cos(2.0*alphaij);


	L = sin(betaij)*cos(alphaij);
	M = sin(betaij)*sin(alphaij);
	N = cos(betaij);


	etapart = (1.0+N*N)*C2a*C2g - 2.0*S2a*S2g*N;
	X = (L*L - M*M) + (eta/3.0)*etapart;
	Y = (3.0*N*N -1.0) + (eta*(1.0-N*N)*C2g);

	etapart = (eta/3.0)*((1+N*N)*S2a*S2g + 2.0*N*C2a*S2g);
	c1s(1,1) = K0 * (2.0*L*M + etapart);			% c -2 case

	etapart = (eta/3.0)*(M*S2g - L*N*C2g);
	c1s(1,2)= K1 * (L*N + etapart);				% c -1 case

	c1s(1,3)= K2 * (sqrt3*X*sin(chi) - Y*cos(chi));		% c  0 case

	etapart = (eta/3.0)*(L* + M*N*C2g);
	c1s(1,4)= K1 * (M*N - etapart);				% c +1 case

	c1s(1,5)= K2 * (sqrt3*X*cos(chi) + Y*sin(chi));		% c +2 case

	% c1s
	% pause

	T1s=zeros(Dim,Dim,5);
	T2s=zeros(Dim,Dim,5);
	kl=1;
	for k=1:n_spins-1
	for l=k+1:n_spins


		if ij==kl

			xi1xi2 = xi1*xi1;

			% Spin tensor for spin i,j change to eigenbasis	
			%{

			T1s(:,:,1)=H_V'*T_D(:,:,1,i,j)*H_V; 
			T1s(:,:,2)=H_V'*T_D(:,:,2,i,j)*H_V; 
			T1s(:,:,3)=H_V'*T_D(:,:,3,i,j)*H_V; 
			T1s(:,:,4)=H_V'*T_D(:,:,4,i,j)*H_V; 
			T1s(:,:,5)=H_V'*T_D(:,:,5,i,j)*H_V; 

			T2s(:,:,1)=H_V'*T_D(:,:,1,i,j)*H_V; 
			T2s(:,:,2)=H_V'*T_D(:,:,2,i,j)*H_V; 
			T2s(:,:,3)=H_V'*T_D(:,:,3,i,j)*H_V; 
			T2s(:,:,4)=H_V'*T_D(:,:,4,i,j)*H_V; 
			T2s(:,:,5)=H_V'*T_D(:,:,5,i,j)*H_V; 

			%}

			T1s(:,:,1)=T_D(:,:,1,i,j); 
			T1s(:,:,2)=T_D(:,:,2,i,j); 
			T1s(:,:,3)=T_D(:,:,3,i,j); 
			T1s(:,:,4)=T_D(:,:,4,i,j); 
			T1s(:,:,5)=T_D(:,:,5,i,j); 

			T2s(:,:,1)=T_D(:,:,1,i,j); 
			T2s(:,:,2)=T_D(:,:,2,i,j); 
			T2s(:,:,3)=T_D(:,:,3,i,j); 
			T2s(:,:,4)=T_D(:,:,4,i,j); 
			T2s(:,:,5)=T_D(:,:,5,i,j); 
			c2s=c1s;


			% xi1xi2
			% c1s
			 % T(:,:,1,i,j)
			%  T(:,:,2,i,j)
			%  T(:,:,3,i,j)
			%  T(:,:,4,i,j)
			%  T(:,:,5,i,j)
			%  pause

			% calculate the Spectral density matrix
			
		for l3=1:5
		mWrf = (l3-2)*Wrf;

			for l1=1:Dim
			for l2=1:Dim	
				J(l1,l2,l3)=0;

				wtr=real(w_omega_rf(1,l1)-w_omega_rf(1,l2))+mWrf;
	%printf("[%d,%d] omega[i]: %0.12f omega[j]: %0.12f w:  %0.12f\n",i,j,real(w_omega_D(1,l1)),real(w_omega_D(1,l2)),wtr);
 
				wtr=wtr*6.283185307;
				for m=1:5
				J_temp=Tau_eff(1,m)/(1+Tau_eff(1,m)*Tau_eff(1,m)*wtr*wtr);
				J(l1,l2,l3)=J(l1,l2,l3)+c1s(1,m)*c2s(1,m)*(J_temp);
	%printf("[%d] tau[i]: %0.12f omega: %0.12f c1[i]: %0.12f c2[i]: %0.12f J: %0.12f\n",m,Tau_eff(1,m),wtr,c1s(1,m),c1s(1,m),J(l1,l2));


				end
				J(l1,l2,l3)=J(l1,l2,l3)/5;
				%printf("J: %0.12f\n",J(l1,l2));
				%pause
				J(l1,l2,l3)=J(l1,l2,l3)*complex(xi1xi2,0);


			end %l1
			end %l2

		end %l3

			% calculate the element of L_DD_RF super operator
			aaa=1;


			Rel=0;
			J2=0;
			J3=0;

			delwaa =0;	% Hetero nuclei secular -omit this for homo
			delwbb =0;	% Hetero nuclei secular -omit this for homo
			for a=1:Dim
			for aa=1:Dim	
				delwaa=w_omega_rf(1,a)-w_omega_rf(1,aa); % Hetero nuclei secular -omit this for homo
				bbb=1;
				for b=1:Dim
				for bb=1:Dim	
				delwbb=w_omega_rf(1,b)-w_omega_rf(1,bb); % Hetero nuclei secular -omit this for homo


				% calculate the relaxa element
				% J2=J(bb,aa);
				% J3=J(a,b);

				g1=0;
				Rel=0;
				k1=1;	


					% printf("J2:%0.10f,%0.10f\n",real(J2),imag(J2));
					% printf("J3:%0.10f,%0.10f\n",real(J3),imag(J3));


					Rel=0;
					for m1=-rank:rank
					Rel=Rel-J(bb,aa,k1)*T1s(a,b,k1)*conj(T2s(aa,bb,k1));

	     				% printf("RII:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(a,b,k1)),imag(T1s(a,b,k1)),real(T2s(aa,bb,k1)),imag(T2s(aa,bb,k1)),real(Rel),imag(Rel));



					Rel=Rel-J(a,b,k1)*T1s(bb,aa,k1)*conj(T2s(b,a,k1));

	    				% printf("RIII:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(bb,aa,k1)),imag(T1s(bb,aa,k1)),real(T2s(b,a,k1)),imag(T2s(b,a,k1)),real(Rel),imag(Rel));


					for g1=1:Dim

					if aa==bb
					Rel=Rel+J(g1,b,k1)*T1s(a,g1,k1)*conj(T2s(b,g1,k1));

				  %  printf("J1:%0.10f,%0.10f\n",real(J(g1,b)),imag(J(g1,b)));
	 			  %  printf("RI:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(a,g1,k1)),imag(T1s(a,g1,k1)),real(T2s(b,g1,k1)),imag(T2s(b,g1,k1)),real(Rel),imag(Rel));



					end

					if a==b
					Rel=Rel+J(bb,g1,k1)*T1s(g1,aa,k1)*conj(T2s(g1,bb,k1));



				  % printf("J1:%0.10f,%0.10f\n",real(J(bb,g1)),imag(J(bb,g1)));
	 			  % printf("RIV:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(g1,aa,k1)),imag(T1s(g1,aa,k1)),real(T2s(g1,bb,k1)),imag(T2s(g1,bb,k1)),real(Rel),imag(Rel));


					end

					end %g1


					k1=k1+1;
					end %m1


				%if Secular_approx ==1 
				%if abs(delwaa-delwbb)<power(10,6) % Secular approx -omit for homo
				%L_DD_RF(aaa,bbb)=L_DD_RF(aaa,bbb)+Rel;
				
				% printf("LR[%d,%d,%d,%d,%d,%d]:%0.6f+%0.6f \n",a,aa,b,bb,aaa,bbb,real(L_DD_RF(aaa,bbb)),imag(L_DD_RF(aaa,bbb)));  					%pause
				%end % fabs(delwaa-delwbb) Secular approx 
				%else
				L_DD_RF(aaa,bbb)=L_DD_RF(aaa,bbb)+Rel;
				%end


				bbb=bbb+1;
				end %bb
				end %b
				% printf("\n");
				aaa=aaa+1;
				%pause
			end %aa
			end %a

			

			%pause
		else % ie ij!=kl

			xi2=xiD(k,l);
			alphakl=Phis(k,l);
			betakl=Thetas(k,l);
			gammakl=0;

			xi1xi2 = xi1*xi2;


			sqrt3 = sqrt(3.0);
			K0 = 0.546274215;		% sqrt[15/(16*pi)]
			K1 = 1.092548431;		% sqrt[15/(4.0*pi)]
			K2 = 0.315391565;		% sqrt[5/(16.0*pi)]

			S2g = sin(2.0*gammakl);
			C2g = cos(2.0*gammakl);
			S2a = sin(2.0*alphakl);
			C2a = cos(2.0*alphakl);


			L = sin(betakl)*cos(alphakl);
			M = sin(betakl)*sin(alphakl);
			N = cos(betakl);


			etapart = (1.0+N*N)*C2a*C2g - 2.0*S2a*S2g*N;
			X = (L*L - M*M) + (eta/3.0)*etapart;
			Y = (3.0*N*N -1.0) + (eta*(1.0-N*N)*C2g);

			etapart = (eta/3.0)*((1+N*N)*S2a*S2g + 2.0*N*C2a*S2g);
			c2s(1,1) = K0 * (2.0*L*M + etapart);			% c -2 case

			etapart = (eta/3.0)*(M*S2g - L*N*C2g);
			c2s(1,2)= K1 * (L*N + etapart);				% c -1 case

			c2s(1,3)= K2 * (sqrt3*X*sin(chi) - Y*cos(chi));		% c  0 case

			etapart = (eta/3.0)*(L* + M*N*C2g);
			c2s(1,4)= K1 * (M*N - etapart);				% c +1 case

			c2s(1,5)= K2 * (sqrt3*X*cos(chi) + Y*sin(chi));		% c +2 case




			% Spin tensor for spin i,j change to eigenbasis	
			%{

			T1s(:,:,1)=H_V'*T_D(:,:,1,i,j)*H_V; 
			T1s(:,:,2)=H_V'*T_D(:,:,2,i,j)*H_V; 
			T1s(:,:,3)=H_V'*T_D(:,:,3,i,j)*H_V; 
			T1s(:,:,4)=H_V'*T_D(:,:,4,i,j)*H_V; 
			T1s(:,:,5)=H_V'*T_D(:,:,5,i,j)*H_V; 

			T2s(:,:,1)=H_V'*T_D(:,:,1,k,l)*H_V; 
			T2s(:,:,2)=H_V'*T_D(:,:,2,k,l)*H_V; 
			T2s(:,:,3)=H_V'*T_D(:,:,3,k,l)*H_V; 
			T2s(:,:,4)=H_V'*T_D(:,:,4,k,l)*H_V; 
			T2s(:,:,5)=H_V'*T_D(:,:,5,k,l)*H_V; 
			%}

			T1s(:,:,1)=T_D(:,:,1,i,j); 
			T1s(:,:,2)=T_D(:,:,2,i,j); 
			T1s(:,:,3)=T_D(:,:,3,i,j); 
			T1s(:,:,4)=T_D(:,:,4,i,j); 
			T1s(:,:,5)=T_D(:,:,5,i,j); 

			T2s(:,:,1)=T_D(:,:,1,k,l); 
			T2s(:,:,2)=T_D(:,:,2,k,l); 
			T2s(:,:,3)=T_D(:,:,3,k,l); 
			T2s(:,:,4)=T_D(:,:,4,k,l); 
			T2s(:,:,5)=T_D(:,:,5,k,l); 
				% xi1xi2
				% c2s
			 	% T(:,:,1,k,l)
				% T(:,:,2,k,l)
				% T(:,:,3,k,l)
				% T(:,:,4,k,l)
				% T(:,:,5,k,l)
			%  pause
			
		for l3=1:5
		mWrf = (l3-2)*Wrf;
		
			% calculate the Spectral density matrix
			for l1=1:Dim
			for l2=1:Dim	
				J(l1,l2,l3)=0;

				wtr=real(w_omega_rf(1,l1)-w_omega_rf(1,l2));
	%printf("[%d,%d] omega[i]: %0.12f omega[j]: %0.12f w:  %0.12f\n",i,j,real(w_omega_D(1,l1)),real(w_omega_D(1,l2)),wtr);
 
				wtr=wtr*6.283185307;
				for m=1:5
				J_temp=Tau_eff(1,m)/(1+Tau_eff(1,m)*Tau_eff(1,m)*wtr*wtr);
				J(l1,l2,l3)=J(l1,l2,l3)+c1s(1,m)*c2s(1,m)*(J_temp);
	%printf("[%d] tau[i]: %0.12f omega: %0.12f c1[i]: %0.12f c2[i]: %0.12f J: %0.12f\n",m,Tau_eff(1,m),wtr,c1s(1,m),c2s(1,m),J(l1,l2));


				end
				J(l1,l2,l3)=J(l1,l2,l3)/5;
				%printf("J: %0.12f\n",J(l1,l2));
				%pause
				J(l1,l2,l3)=J(l1,l2,l3)*complex(xi1xi2,0);


			end %l1
			end %l2
		end %l3


			% calculate the element of L_DD_RF super operator
			aaa=1;


			Rel=0;
			J2=0;
			J3=0;
			delwaa=0;
			delwbb=0;

			for a=1:Dim
			for aa=1:Dim
					delwaa = w_omega_rf(1,a)-w_omega_rf(aa);			
				bbb=1;
				for b=1:Dim
				for bb=1:Dim	


					delwbb = w_omega_rf(1,b)-w_omega_rf(bb);
				% calculate the relaxa element
				%J2=J(bb,aa);
				%J3=J(a,b);

				g1=0;
				Rel=0;
				k1=1;	


					% printf("J2:%0.10f,%0.10f\n",real(J2),imag(J2));
					% printf("J3:%0.10f,%0.10f\n",real(J3),imag(J3));


					Rel=0;
					for m1=-rank:rank
					Rel=Rel-J(bb,aa,k1)*T1s(a,b,k1)*conj(T2s(aa,bb,k1));

	     				% printf("RII:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(a,b,k1)),imag(T1s(a,b,k1)),real(T2s(aa,bb,k1)),imag(T2s(aa,bb,k1)),real(Rel),imag(Rel));



					Rel=Rel-J(a,b,k1)*T1s(bb,aa,k1)*conj(T2s(b,a,k1));

	    				% printf("RIII:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(bb,aa,k1)),imag(T1s(bb,aa,k1)),real(T2s(b,a,k1)),imag(T2s(b,a,k1)),real(Rel),imag(Rel));


					for g1=1:Dim

					if aa==bb
					Rel=Rel+J(g1,b,k1)*T1s(a,g1,k1)*conj(T2s(b,g1,k1));

				  %  printf("J1:%0.10f,%0.10f\n",real(J(g1,b)),imag(J(g1,b)));
	 			  %  printf("RI:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(a,g1,k1)),imag(T1s(a,g1,k1)),real(T2s(b,g1,k1)),imag(T2s(b,g1,k1)),real(Rel),imag(Rel));



					end

					if a==b
					Rel=Rel+J(bb,g1,k1)*T1s(g1,aa,k1)*conj(T2s(g1,bb,k1));



				  % printf("J1:%0.10f,%0.10f\n",real(J(bb,g1)),imag(J(bb,g1)));
	 			  % printf("RIV:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(g1,aa,k1)),imag(T1s(g1,aa,k1)),real(T2s(g1,bb,k1)),imag(T2s(g1,bb,k1)),real(Rel),imag(Rel));


					end

					end %g1


					k1=k1+1;
					end %m1

				%if Secular_approx ==1 
				%if abs(delwaa-delwbb)<power(10,6) % Secular approx -omit for homo
				%L_DD_RF(aaa,bbb)=L_DD_RF(aaa,bbb)+Rel;
				
				% printf("LR[%d,%d,%d,%d,%d,%d]:%0.6f+%0.6f \n",a,aa,b,bb,aaa,bbb,real(L_DD_RF(aaa,bbb)),imag(L_DD_RF(aaa,bbb)));  					%pause
				%end % fabs(delwaa-delwbb) Secular approx 
				%else
				L_DD_RF(aaa,bbb)=L_DD_RF(aaa,bbb)+Rel;
				%end

				bbb=bbb+1;
				end %bb
				end %b
				% printf("\n");
				aaa=aaa+1;
				%pause
			end %aa
			end %a

			

			%pause


		end %if ij==kl



	%printf("ij[%d %d  %d  %d]: %d \n",i,j,k,l,ij);
	%printf("kl[%d %d  %d  %d]: %d \n",i,j,k,l,kl);


	kl=kl+1;
	end %l
	end %k

ij=ij+1;
end %j
end %i






%{
%-----------------%
% CSA mechanism	  %
%-----------------%

% Calculate diffusion constants (See DD ): diff_x,y,z;Tau_eff
% Calculate Chi constant (See DD ): chi

% setup HC for CSA in Lab frame  (used by DD-CSA)

HC=-(v_Hz_offseted(1,1)+MHZ_H*power(10,6)*(Gamma_1(1,1)/Gamma_H))*HzNi - (v_Hz_offseted(1,2)+MHZ_H*power(10,6)*Gamma_1(1,2)/Gamma_H)*HiNz;

HC=HC+T_Jham;
%HC=H_V'*HC*H_V;





% w_omega for CSA
w_omega_CS=zeros(1,Dim);
for i=1:Dim
w_omega_CS(1,i)=real(HC(i,i));
%w_omega_CS(1,i)=real(HC_D(i,i));
%printf("w_omega_CS[%d]:%f (HZ:%f  Ho:%f)\n",i,w_omega_CS(1,i),real(HZ(i,i)),real(H_D(i,i)));
end



% CSA interaction constants between all spin pairs by calculating inter nuclear distances (used by DD-CSA)

xiC=zeros(n_spins,n_spins);
a=1.941625913;
%K = -2.0*a*mu0d4pi*hbar;
K = 1.941625913;  %sqrt(6*pi/5)
omiomi=0;
delzz=0;
xii=0;
for i=1:n_spins
	omi=MHZ(1,i);
	delzz=2/3*v_aniso(1,i);
	xiC(i,i)=K*pi2*omi*delzz;
   %printf("[%d] K:%f omi:%f delzz:%f xiC:%f\n",i,K,omi,delzz,xiC(i,i)); pause

end





% Generate the spin tensor in AAS for CSA (used by DD-CSA)

T1_CS_20=zeros(Dim,Dim);
T1_CS_2_1=zeros(Dim,Dim);
T1_CS_2_2=zeros(Dim,Dim);
T1_CS_21=zeros(Dim,Dim);
T1_CS_22=zeros(Dim,Dim);

T2_CS_20=zeros(Dim,Dim);
T2_CS_2_1=zeros(Dim,Dim);
T2_CS_2_2=zeros(Dim,Dim);
T2_CS_21=zeros(Dim,Dim);
T2_CS_22=zeros(Dim,Dim);

B0=[0 0 1]; % along z axis
Bx=B0(1,1);
By=B0(1,2);
Bz=B0(1,3);
Bp=complex(Bx,By);
Bm=complex(Bx,-By);

T1_CS_20=1/sqrt(6)*(3*HzNi.*Bz-(HxNi*Bx+HyNi*By+HzNi*Bz));
T1_CS_21=-0.5*(HzNi.*Bp+HpNi.*Bz);
T1_CS_2_1=0.5*(HzNi.*Bm+HmNi.*Bz);
T1_CS_22=0.5*(HpNi.*Bp);
T1_CS_2_2=0.5*(HmNi.*Bm);


T2_CS_20=1/sqrt(6)*(3*HiNz.*Bz-(HiNx*Bx+HiNy*By+HiNz*Bz));
T2_CS_21=-0.5*(HiNz.*Bp+HiNp.*Bz);
T2_CS_2_1=0.5*(HiNz.*Bm+HiNm.*Bz);
T2_CS_22=0.5*(HiNp.*Bp);
T2_CS_2_2=0.5*(HiNm.*Bm);

% Put all the spin tensors in a 3d array with last two indices as (i,j): total 4D array (used by DD-CSA)
T_CS=zeros(Dim,Dim,5,n_spins);

for i=1:2
	if i==1
		T_CS(:,:,1,i)=T1_CS_2_2(:,:);
		T_CS(:,:,2,i)=T1_CS_2_1(:,:);
		T_CS(:,:,3,i)=T1_CS_20(:,:);
		T_CS(:,:,4,i)=T1_CS_21(:,:);
		T_CS(:,:,5,i)=T1_CS_22(:,:);
	elseif i==2
		T_CS(:,:,1,i)=T2_CS_2_2(:,:);
		T_CS(:,:,2,i)=T2_CS_2_1(:,:);
		T_CS(:,:,3,i)=T2_CS_20(:,:);
		T_CS(:,:,4,i)=T2_CS_21(:,:);
		T_CS(:,:,5,i)=T2_CS_22(:,:);
	end %i==j
end %i


% Refield matrix element

xi1=0;
xi2=0;
xi1xi2=0;

% Dipolar variables
alphaij=0;
betaij=0;
alphakl=0;
betakl=0;
ij=1;
kl=0;
w=0;
w0=0;
w1=0;
w2=0;
wi=0;
wj=0;
eta=0;

% CSA variables
alphai=0;
alphaj=0;
betai=0;
betaj=0;
gammai=0;
gammaj=0;
eta=0;
sqrt3 =0;
K0 = 0;
K1 = 0;
K2 = 0;
S2g = 0;
C2g = 0;
S2a = 0;
C2a = 0;
L = 0;
M = 0;
N = 0;
etapart = 0;
X = 0;
Y = 0;

% Common to both
rank=2;
J=zeros(Dim,Dim);
J_temp=0;
EA1=zeros(1,3);
EA2=zeros(1,3);
c1s=zeros(1,5);
c2s=zeros(1,5);
sqrt3 = sqrt(3.0);


%for i=1:n_spins-1	% Dipolar
%for j=i+1:n_spins	% Dipolar

for i=1:n_spins		% CSA
for j=1:n_spins		% CSA
	%{
	xi1=xiD(i,j);	% Dipolar
	alphaij=Phis(i,j);
	betaij=Thetas(i,j);
	gammaij=0;
	%}


	xi1=real(xiC(i,i)); % CSA
	alphai=EA(i,1);
	betai=EA(i,2);
	gammai=EA(i,3);

	% spectral density and spatial spherical harmonics coefficient combined for 5 elements of spatial operator :: five spectral density function coefficients  for isotropic rigid diffusion
	% EA   : Euler angles (radians)
	% chi  : Ratio of rotational correlation times
	% eta  : Tensor value

	% Same for both mechanism just change the angles from EA for CSA, EA calculated from coord for dipolar

	alphaij=alphai;
	betaij=betai;
	gammaij=gammai;

	K0 = 0.546274215;		% sqrt[15/(16*pi)]
	K1 = 1.092548431;		% sqrt[15/(4.0*pi)]
	K2 = 0.315391565;		% sqrt[5/(16.0*pi)]

	S2g = sin(2.0*gammaij);
	C2g = cos(2.0*gammaij);
	S2a = sin(2.0*alphaij);
	C2a = cos(2.0*alphaij);


	L = sin(betaij)*cos(alphaij);
	M = sin(betaij)*sin(alphaij);
	N = cos(betaij);


	etapart = (1.0+N*N)*C2a*C2g - 2.0*S2a*S2g*N;
	X = (L*L - M*M) + (eta/3.0)*etapart;
	Y = (3.0*N*N -1.0) + (eta*(1.0-N*N)*C2g);

	etapart = (eta/3.0)*((1+N*N)*S2a*S2g + 2.0*N*C2a*S2g);
	c1s(1,1) = K0 * (2.0*L*M + etapart);			% c -2 case

	etapart = (eta/3.0)*(M*S2g - L*N*C2g);
	c1s(1,2)= K1 * (L*N + etapart);				% c -1 case

	c1s(1,3)= K2 * (sqrt3*X*sin(chi) - Y*cos(chi));		% c  0 case

	etapart = (eta/3.0)*(L* + M*N*C2g);
	c1s(1,4)= K1 * (M*N - etapart);				% c +1 case

	c1s(1,5)= K2 * (sqrt3*X*cos(chi) + Y*sin(chi));		% c +2 case

	% c1s
	% pause

	T1s=zeros(Dim,Dim,5);
	T2s=zeros(Dim,Dim,5);
	%kl=1; % Dipolar

	% for k=1:n_spins-1 % dipolar 3rd
	% for l=k+1:n_spins % dipolar 4th


		% if ij==kl	% dipolar
		if i==j	% CSA
			xi1xi2 = xi1*xi1;

			% Spin tensor for spin i,j change to eigenbasis	

			%{
			T1s(:,:,1)=H_V'*T_CS(:,:,1,i)*H_V; 
			T1s(:,:,2)=H_V'*T_CS(:,:,2,i)*H_V; 
			T1s(:,:,3)=H_V'*T_CS(:,:,3,i)*H_V; 
			T1s(:,:,4)=H_V'*T_CS(:,:,4,i)*H_V; 
			T1s(:,:,5)=H_V'*T_CS(:,:,5,i)*H_V; 

			T2s(:,:,1)=H_V'*T_CS(:,:,1,j)*H_V; 
			T2s(:,:,2)=H_V'*T_CS(:,:,2,j)*H_V; 
			T2s(:,:,3)=H_V'*T_CS(:,:,3,j)*H_V; 
			T2s(:,:,4)=H_V'*T_CS(:,:,4,j)*H_V; 
			T2s(:,:,5)=H_V'*T_CS(:,:,5,j)*H_V; 
			%}
			
			T1s(:,:,1)=T_CS(:,:,1,i); 
			T1s(:,:,2)=T_CS(:,:,2,i); 
			T1s(:,:,3)=T_CS(:,:,3,i); 
			T1s(:,:,4)=T_CS(:,:,4,i); 
			T1s(:,:,5)=T_CS(:,:,5,i); 

			T2s(:,:,1)=T_CS(:,:,1,j); 
			T2s(:,:,2)=T_CS(:,:,2,j); 
			T2s(:,:,3)=T_CS(:,:,3,j); 
			T2s(:,:,4)=T_CS(:,:,4,j); 
			T2s(:,:,5)=T_CS(:,:,5,j); 


			c2s=c1s;
			% xi1xi2
			% c1s
			 % T(:,:,1,i,j)
			%  T(:,:,2,i,j)
			%  T(:,:,3,i,j)
			%  T(:,:,4,i,j)
			%  T(:,:,5,i,j)
			%  pause

			% calculate the Spectral density matrix same for any mechanism depends only on CS and tau
			for l1=1:Dim
			for l2=1:Dim	
				J(l1,l2)=0;

				wtr=real(w_omega_CS(1,l1)-w_omega_CS(1,l2));
	%printf("[%d,%d] omega[i]: %0.12f omega[j]: %0.12f w:  %0.12f\n",i,j,real(w_omega_CS(1,l1)),real(w_omega_CS(1,l2)),wtr);
 
				wtr=wtr*6.283185307;
				for m=1:5
				J_temp=Tau_eff(1,m)/(1+Tau_eff(1,m)*Tau_eff(1,m)*wtr*wtr);
				J(l1,l2)=J(l1,l2)+c1s(1,m)*c2s(1,m)*(J_temp);
	%printf("[%d] tau[i]: %0.12f omega: %0.12f c1[i]: %0.12f c2[i]: %0.12f J: %0.12f\n",m,Tau_eff(1,m),wtr,c1s(1,m),c1s(1,m),J(l1,l2));


				end
				J(l1,l2)=J(l1,l2)/5;
				%printf("J: %0.12f\n",J(l1,l2));
				%pause
				J(l1,l2)=J(l1,l2)*complex(xi1xi2,0);


			end %l1
			end %l2

			%printf("\n\n\nLoop[%d][%d]:\n",i,j); pause
			%xi1xi2
			%T1s
			%T2s
			%c1s
			%c1s
			%printf("\nJ\n");
			%for l1=1:Dim
			%for l2=1:Dim	
			%	printf("%0.10f ",J(l1,l2));
			%end
			%printf("\n");
			%end
			%pause 			
			%printf("\n\n\n");


			% calculate the element of L_DD_RF super operator
			aaa=1;


			Rel=0;
			J2=0;
			J3=0;


			for a=1:Dim
			for aa=1:Dim	
				delwaa=w_omega_CS(1,a)-w_omega_CS(1,aa); % Hetero nuclei secular -omit this for homo

				bbb=1;
				for b=1:Dim
				for bb=1:Dim	
				delwbb=w_omega_CS(1,b)-w_omega_CS(1,bb); % Hetero nuclei secular -omit this for homo

				% calculate the relaxa element
				J2=J(bb,aa);
				J3=J(a,b);

				g1=0;
				k1=1;	


					% printf("J2:%0.10f,%0.10f\n",real(J2),imag(J2));
					% printf("J3:%0.10f,%0.10f\n",real(J3),imag(J3));


					Rel=0;
					for m1=-rank:rank
					Rel=Rel-J2*T1s(a,b,k1)*T2s(aa,bb,k1);

	     				% printf("RII:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(a,b,k1)),imag(T1s(a,b,k1)),real(T2s(aa,bb,k1)),imag(T2s(aa,bb,k1)),real(Rel),imag(Rel));



					Rel=Rel-J3*T1s(bb,aa,k1)*T2s(b,a,k1);

	    				% printf("RIII:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(bb,aa,k1)),imag(T1s(bb,aa,k1)),real(T2s(b,a,k1)),imag(T2s(b,a,k1)),real(Rel),imag(Rel));


					for g1=1:Dim

					if aa==bb
					Rel=Rel+J(g1,b)*T1s(a,g1,k1)*T2s(b,g1,k1);

				  % printf("J1:%0.10f,%0.10f\n",real(J(g1,b)),imag(J(g1,b)));
	 			  % printf("RI:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(a,g1,k1)),imag(T1s(a,g1,k1)),real(T2s(b,g1,k1)),imag(T2s(b,g1,k1)),real(Rel),imag(Rel));



					end

					if a==b
					Rel=Rel+J(bb,g1)*T1s(g1,aa,k1)*T2s(g1,bb,k1);



				  % printf("J1:%0.10f,%0.10f\n",real(J(bb,g1)),imag(J(bb,g1)));
	 			 %  printf("RIV:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(g1,aa,k1)),imag(T1s(g1,aa,k1)),real(T2s(g1,bb,k1)),imag(T2s(g1,bb,k1)),real(Rel),imag(Rel));


					end

					end %g1


					k1=k1+1;
					end %m1


				if Secular_approx ==1 
				if abs(delwaa-delwbb)<power(10,6) % Secular approx -omit for homo
				L_CC(aaa,bbb)=L_CC(aaa,bbb)+Rel;
				end
				else
				L_CC(aaa,bbb)=L_CC(aaa,bbb)+Rel;
				% printf("LR[%d,%d,%d,%d,%d,%d]:%0.8f+%0.8f \n",a,aa,b,bb,aaa,bbb,real(L_CC(aaa,bbb)),imag(L_CC(aaa,bbb)));  				pause
			%printf("LR[%d,%d,%d,%d]:%0.8f+%0.8f Cum:%0.8f+%0.8f \n",a,aa,b,bb,real(Rel),imag(Rel),real(L_CC(aaa,bbb)),imag(L_CC(aaa,bbb)));  		

				end % 	Secular_approx	
				bbb=bbb+1;
				end %bb
				end %b
				% printf("\n");				pause
				aaa=aaa+1;
			%printf("\n");
			end %aa
			%printf("\n");
			end %a
			%printf("\n");
			%pause
			%L_CC
			% pause
			else 	% else ij~=kl

			%{
			xi2=xiD(k,l);	%Dipolar
			alphakl=Phis(k,l);
			betakl=Thetas(k,l);
			gammakl=0;
			%}
			% Same for both mechanism just change the angles and index to j in EA for CSA, EA calculated from coord for dipolar

			xi2=real(xiC(j,j)); %CSA

			xi1xi2 = xi1*xi2;

			alphakl= EA(j,1);
			betakl= EA(j,2);
			gammakl= EA(j,3);


			K0 = 0.546274215;		% sqrt[15/(16*pi)]
			K1 = 1.092548431;		% sqrt[15/(4.0*pi)]
			K2 = 0.315391565;		% sqrt[5/(16.0*pi)]

			S2g = sin(2.0*gammakl);
			C2g = cos(2.0*gammakl);
			S2a = sin(2.0*alphakl);
			C2a = cos(2.0*alphakl);


			L = sin(betakl)*cos(alphakl);
			M = sin(betakl)*sin(alphakl);
			N = cos(betakl);


			etapart = (1.0+N*N)*C2a*C2g - 2.0*S2a*S2g*N;
			X = (L*L - M*M) + (eta/3.0)*etapart;
			Y = (3.0*N*N -1.0) + (eta*(1.0-N*N)*C2g);

			etapart = (eta/3.0)*((1+N*N)*S2a*S2g + 2.0*N*C2a*S2g);
			c2s(1,1) = K0 * (2.0*L*M + etapart);			% c -2 case

			etapart = (eta/3.0)*(M*S2g - L*N*C2g);
			c2s(1,2)= K1 * (L*N + etapart);				% c -1 case

			c2s(1,3)= K2 * (sqrt3*X*sin(chi) - Y*cos(chi));		% c  0 case

			etapart = (eta/3.0)*(L* + M*N*C2g);
			c2s(1,4)= K1 * (M*N - etapart);				% c +1 case

			c2s(1,5)= K2 * (sqrt3*X*cos(chi) + Y*sin(chi));		% c +2 case




			% Spin tensor for spin i,j change to eigenbasis	for Dipolar

			%{
			T1s(:,:,1)=H_V'*T_CS(:,:,1,i)*H_V; 
			T1s(:,:,2)=H_V'*T_CS(:,:,2,i)*H_V; 
			T1s(:,:,3)=H_V'*T_CS(:,:,3,i)*H_V; 
			T1s(:,:,4)=H_V'*T_CS(:,:,4,i)*H_V; 
			T1s(:,:,5)=H_V'*T_CS(:,:,5,i)*H_V; 

			T2s(:,:,1)=H_V'*T_CS(:,:,1,j)*H_V; 
			T2s(:,:,2)=H_V'*T_CS(:,:,2,j)*H_V; 
			T2s(:,:,3)=H_V'*T_CS(:,:,3,j)*H_V; 
			T2s(:,:,4)=H_V'*T_CS(:,:,4,j)*H_V; 
			T2s(:,:,5)=H_V'*T_CS(:,:,5,j)*H_V; 
			%}

			T1s(:,:,1)=T_CS(:,:,1,i); 
			T1s(:,:,2)=T_CS(:,:,2,i); 
			T1s(:,:,3)=T_CS(:,:,3,i); 
			T1s(:,:,4)=T_CS(:,:,4,i);
			T1s(:,:,5)=T_CS(:,:,5,i); 

			T2s(:,:,1)=T_CS(:,:,1,j); 
			T2s(:,:,2)=T_CS(:,:,2,j); 
			T2s(:,:,3)=T_CS(:,:,3,j); 
			T2s(:,:,4)=T_CS(:,:,4,j); 
			T2s(:,:,5)=T_CS(:,:,5,j); 
				% xi1xi2
				% c2s
			 	% T(:,:,1,k,l)
				% T(:,:,2,k,l)
				% T(:,:,3,k,l)
				% T(:,:,4,k,l)
				% T(:,:,5,k,l)
			%  pause

			% calculate the Spectral density matrix
			for l1=1:Dim
			for l2=1:Dim	
				J(l1,l2)=0;

				wtr=real(w_omega_CS(1,l1)-w_omega_CS(1,l2));
	%printf("[%d,%d] omega[i]: %0.12f omega[j]: %0.12f w:  %0.12f\n",i,j,real(w_omega_CS(1,l1)),real(w_omega_CS(1,l2)),wtr);
 
				wtr=wtr*6.283185307;
				for m=1:5
				J_temp=Tau_eff(1,m)/(1+Tau_eff(1,m)*Tau_eff(1,m)*wtr*wtr);
				J(l1,l2)=J(l1,l2)+c1s(1,m)*c2s(1,m)*(J_temp);
	%printf("[%d] tau[i]: %0.12f omega: %0.12f c1[i]: %0.12f c2[i]: %0.12f J: %0.12f\n",m,Tau_eff(1,m),wtr,c1s(1,m),c2s(1,m),J(l1,l2));


				end
				J(l1,l2)=J(l1,l2)/5;
				%printf("J: %0.12f\n",J(l1,l2));
				%pause
				J(l1,l2)=J(l1,l2)*complex(xi1xi2,0);


			end %l1
			end %l2

			%  printf("\n\n\nLoop[%d][%d]:\n",i,j); pause
			%  xi1xi2
			%  T1s
			%  T2s
			%  c1s
			%  c2s
			%  J
			%  pause

			% calculate the element of L_DD_RF super operator
			aaa=1;


			Rel=0;
			J2=0;
			J3=0;


			for a=1:Dim
			for aa=1:Dim	

				delwaa=w_omega_CS(1,a)-w_omega_CS(1,aa); % Hetero nuclei secular -omit this for homo

				bbb=1;
				for b=1:Dim
				for bb=1:Dim	
				delwbb=w_omega_CS(1,b)-w_omega_CS(1,bb); % Hetero nuclei secular -omit this for homo

				% calculate the relaxa element
				J2=J(bb,aa);
				J3=J(a,b);

				g1=0;
				Rel=0;
				k1=1;	


					% printf("J2:%0.10f,%0.10f\n",real(J2),imag(J2));
					% printf("J3:%0.10f,%0.10f\n",real(J3),imag(J3));


					Rel=0;
					for m1=-rank:rank
					Rel=Rel-J2*T1s(a,b,k1)*T2s(aa,bb,k1);

	     				% printf("RII:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(a,b,k1)),imag(T1s(a,b,k1)),real(T2s(aa,bb,k1)),imag(T2s(aa,bb,k1)),real(Rel),imag(Rel));



					Rel=Rel-J3*T1s(bb,aa,k1)*T2s(b,a,k1);

	    				% printf("RIII:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(bb,aa,k1)),imag(T1s(bb,aa,k1)),real(T2s(b,a,k1)),imag(T2s(b,a,k1)),real(Rel),imag(Rel));


					for g1=1:Dim

					if aa==bb
					Rel=Rel+J(g1,b)*T1s(a,g1,k1)*T2s(b,g1,k1);

				  %  printf("J1:%0.10f,%0.10f\n",real(J(g1,b)),imag(J(g1,b)));
	 			  %  printf("RI:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(a,g1,k1)),imag(T1s(a,g1,k1)),real(T2s(b,g1,k1)),imag(T2s(b,g1,k1)),real(Rel),imag(Rel));



					end

					if a==b
					Rel=Rel+J(bb,g1)*T1s(g1,aa,k1)*T2s(g1,bb,k1);



				  % printf("J1:%0.10f,%0.10f\n",real(J(bb,g1)),imag(J(bb,g1)));
	 			  % printf("RIV:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(g1,aa,k1)),imag(T1s(g1,aa,k1)),real(T2s(g1,bb,k1)),imag(T2s(g1,bb,k1)),real(Rel),imag(Rel));


					end

					end %g1


					k1=k1+1;
					end %m1

				if Secular_approx ==1 
				if abs(delwaa-delwbb)<power(10,6) % Secular approx -omit for homo
				L_CC(aaa,bbb)=L_CC(aaa,bbb)+Rel;
				end
				else
				L_CC(aaa,bbb)=L_CC(aaa,bbb)+Rel;
				% printf("LR[%d,%d,%d,%d,%d,%d]:%0.8f+%0.8f \n",a,aa,b,bb,aaa,bbb,real(L_CC(aaa,bbb)),imag(L_CC(aaa,bbb)));  				pause
			%printf("LR[%d,%d,%d,%d]:%0.8f+%0.8f Cum:%0.8f+%0.8f \n",a,aa,b,bb,real(Rel),imag(Rel),real(L_CC(aaa,bbb)),imag(L_CC(aaa,bbb)));  		

				end % 	Secular_approx	

				bbb=bbb+1;
				end %bb
				end %b
				% printf("\n");
				aaa=aaa+1;
				%pause
			%  printf("\n");
			end %aa
			%  printf("\n");
			end %a

			%  printf("\n");
			%pause
			%L_CC
			%pause
			%pause


		end %if ij==kl



	%printf("ij[%d %d  %d  %d]: %d \n",i,j,k,l,ij);
	%printf("kl[%d %d  %d  %d]: %d \n",i,j,k,l,kl);


	% kl=kl+1; % Dipolar

	% end %l 3rd loop for dipolar
	% end %k 4th loop for dipolar

% ij=ij+1; % Dipolar
end %j
end %i


%-----------------------%
% Dipolar-CSA mechanism	%
%-----------------------%

% Calculate diffusion constants (See DD): diff_x,y,z;Tau_eff

% Calculate Chi constant (See DD): chi


% setup HC for CSA in Lab frame  (See CSA)

% CSA interaction constants between all spin pairs by calculating inter nuclear distances (See CSA)

% Dipolar interaction constants between all spin pairs by calculating inter nuclear distances

% Note: XiD required is precalculated in DD mechanism using it here

% CSA interaction constants between all spin pairs by calculating inter nuclear distances

% Generate the spin tensor in AAS for dipolar

% Spin tensor pre-calc in DD mechanism is required here  {T12_D_*}

% Put all the spin tensors in a 3d array with last two indices as (i,j): total 5D array {T_D}

% Generate the spin tensor in AAS for CSA (See - CSA)

% Put all the spin tensors in a 3d array with last two indices as (i,j): total 4D array (See CSA)

% Refield matrix element

xi1=0;
xi2=0;
xi1xi2=0;

% Dipolar variables
alphaij=0;
betaij=0;
alphakl=0;
betakl=0;
ij=1;
kl=0;
w=0;
w0=0;
w1=0;
w2=0;
wi=0;
wj=0;
eta=0;

% CSA variables
alphai=0;
alphaj=0;
betai=0;
betaj=0;
gammai=0;
gammaj=0;
eta=0;
sqrt3 =0;
K0 = 0;
K1 = 0;
K2 = 0;
S2g = 0;
C2g = 0;
S2a = 0;
C2a = 0;
L = 0;
M = 0;
N = 0;
etapart = 0;
X = 0;
Y = 0;

% Common to both
rank=2;
J=zeros(Dim,Dim);
J_temp=0;
EA1=zeros(1,3);
EA2=zeros(1,3);
c1s=zeros(1,5);
c2s=zeros(1,5);
sqrt3 = sqrt(3.0);

% First cross interaction between dipolarloop ij vs CSA loop k

for i=1:n_spins-1	% Dipolar
for j=i+1:n_spins	% Dipolar

	xi1=xiD(i,j);	% Dipolar
	alphaij=Phis(i,j);
	betaij=Thetas(i,j);
	gammaij=0;

	K0 = 0.546274215;		% sqrt[15/(16*pi)]
	K1 = 1.092548431;		% sqrt[15/(4.0*pi)]
	K2 = 0.315391565;		% sqrt[5/(16.0*pi)]

	S2g = sin(2.0*gammaij);
	C2g = cos(2.0*gammaij);
	S2a = sin(2.0*alphaij);
	C2a = cos(2.0*alphaij);


	L = sin(betaij)*cos(alphaij);
	M = sin(betaij)*sin(alphaij);
	N = cos(betaij);


	etapart = (1.0+N*N)*C2a*C2g - 2.0*S2a*S2g*N;
	X = (L*L - M*M) + (eta/3.0)*etapart;
	Y = (3.0*N*N -1.0) + (eta*(1.0-N*N)*C2g);

	etapart = (eta/3.0)*((1+N*N)*S2a*S2g + 2.0*N*C2a*S2g);
	c1s(1,1) = K0 * (2.0*L*M + etapart);			% c -2 case

	etapart = (eta/3.0)*(M*S2g - L*N*C2g);
	c1s(1,2)= K1 * (L*N + etapart);				% c -1 case

	c1s(1,3)= K2 * (sqrt3*X*sin(chi) - Y*cos(chi));		% c  0 case

	etapart = (eta/3.0)*(L* + M*N*C2g);
	c1s(1,4)= K1 * (M*N - etapart);				% c +1 case

	c1s(1,5)= K2 * (sqrt3*X*cos(chi) + Y*sin(chi));		% c +2 case

	% c1s
	% pause


			T1s=zeros(Dim,Dim,5);


			%{
			T1s(:,:,1)=H_V'*T_D(:,:,1,i,j)*H_V; 
			T1s(:,:,2)=H_V'*T_D(:,:,2,i,j)*H_V; 
			T1s(:,:,3)=H_V'*T_D(:,:,3,i,j)*H_V; 
			T1s(:,:,4)=H_V'*T_D(:,:,4,i,j)*H_V; 
			T1s(:,:,5)=H_V'*T_D(:,:,5,i,j)*H_V; 
			%}

			T1s(:,:,1)=T_D(:,:,1,i,j); 
			T1s(:,:,2)=T_D(:,:,2,i,j); 
			T1s(:,:,3)=T_D(:,:,3,i,j); 
			T1s(:,:,4)=T_D(:,:,4,i,j); 
			T1s(:,:,5)=T_D(:,:,5,i,j); 

			%{
			T1s(:,:,1)=T_D(:,:,1,i,j); 
			T1s(:,:,2)=T_D(:,:,2,i,j); 
			T1s(:,:,3)=T_D(:,:,3,i,j); 
			T1s(:,:,4)=T_D(:,:,4,i,j); 
			T1s(:,:,5)=T_D(:,:,5,i,j); 
			%}

		for k=1:n_spins	% Dipolar

			xi2=real(xiC(k,k)); % CSA
			alphai=EA(k,1);
			betai=EA(k,2);
			gammai=EA(k,3);

			alphaij=alphai;
			betaij=betai;
			gammaij=gammai;

			K0 = 0.546274215;		% sqrt[15/(16*pi)]
			K1 = 1.092548431;		% sqrt[15/(4.0*pi)]
			K2 = 0.315391565;		% sqrt[5/(16.0*pi)]

			S2g = sin(2.0*gammaij);
			C2g = cos(2.0*gammaij);
			S2a = sin(2.0*alphaij);
			C2a = cos(2.0*alphaij);


			L = sin(betaij)*cos(alphaij);
			M = sin(betaij)*sin(alphaij);
			N = cos(betaij);


			etapart = (1.0+N*N)*C2a*C2g - 2.0*S2a*S2g*N;
			X = (L*L - M*M) + (eta/3.0)*etapart;
			Y = (3.0*N*N -1.0) + (eta*(1.0-N*N)*C2g);

			etapart = (eta/3.0)*((1+N*N)*S2a*S2g + 2.0*N*C2a*S2g);
			c2s(1,1) = K0 * (2.0*L*M + etapart);			% c -2 case

			etapart = (eta/3.0)*(M*S2g - L*N*C2g);
			c2s(1,2)= K1 * (L*N + etapart);				% c -1 case

			c2s(1,3)= K2 * (sqrt3*X*sin(chi) - Y*cos(chi));		% c  0 case

			etapart = (eta/3.0)*(L* + M*N*C2g);
			c2s(1,4)= K1 * (M*N - etapart);				% c +1 case

			c2s(1,5)= K2 * (sqrt3*X*cos(chi) + Y*sin(chi));		% c +2 case


			T2s=zeros(Dim,Dim,5);
			%{
			T2s(:,:,1)=H_V'*T_CS(:,:,1,k)*H_V; 
			T2s(:,:,2)=H_V'*T_CS(:,:,2,k)*H_V; 
			T2s(:,:,3)=H_V'*T_CS(:,:,3,k)*H_V; 
			T2s(:,:,4)=H_V'*T_CS(:,:,4,k)*H_V; 
			T2s(:,:,5)=H_V'*T_CS(:,:,5,k)*H_V; 
			%}
			
			T2s(:,:,1)=T_CS(:,:,1,k); 
			T2s(:,:,2)=T_CS(:,:,2,k); 
			T2s(:,:,3)=T_CS(:,:,3,k); 
			T2s(:,:,4)=T_CS(:,:,4,k); 
			T2s(:,:,5)=T_CS(:,:,5,k); 

			%{
			T2s(:,:,1)=T_CS(:,:,1,k); 
			T2s(:,:,2)=T_CS(:,:,2,k); 
			T2s(:,:,3)=T_CS(:,:,3,k); 
			T2s(:,:,4)=T_CS(:,:,4,k); 
			T2s(:,:,5)=T_CS(:,:,5,k); 
			%}

			
			% calculate the Spectral density matrix same for any mechanism depends only on CS and tau
			xi1xi2=xi1*xi2;
			for l1=1:Dim
			for l2=1:Dim	
				J(l1,l2)=0;

				wtr=real(w_omega_D(1,l1)-w_omega_D(1,l2));
	%printf("[%d,%d] omega[i]: %0.12f omega[j]: %0.12f w:  %0.12f\n",i,j,real(w_omega(1,l1)),real(w_omega(1,l2)),wtr);
 
				wtr=wtr*6.283185307;
				for m=1:5
				J_temp=Tau_eff(1,m)/(1+Tau_eff(1,m)*Tau_eff(1,m)*wtr*wtr);
				J(l1,l2)=J(l1,l2)+c1s(1,m)*c2s(1,m)*(J_temp);
	%printf("[%d] tau[i]: %0.12f omega: %0.12f c1[i]: %0.12f c2[i]: %0.12f J: %0.12f\n",m,Tau_eff(1,m),wtr,c1s(1,m),c1s(1,m),J(l1,l2));


				end
				J(l1,l2)=J(l1,l2)/5;
				%printf("J: %0.12f\n",J(l1,l2));
				%pause
				J(l1,l2)=J(l1,l2)*complex(xi1xi2,0);


			end %l1
			end %l2


			%{
			printf("\n\n\nLoop[%d][%d][%d]:\n",i,j,k); pause
			xi1xi2
			T1s
			T2s
			c1s
			c2s
			printf("\nJ\n");
			for l1=1:Dim
			for l2=1:Dim	
				printf("%0.10f ",J(l1,l2));
			end
			printf("\n");
			end
			w_omega_D
			 pause 			
			% printf("\n\n\n");
			%}

			% calculate the element of L_DC super operator
			aaa=1;


			Rel=0;
			J2=0;
			J3=0;


			for a=1:Dim
			for aa=1:Dim	
				delwaa=w_omega_D(1,a)-w_omega_D(1,aa); % Hetero nuclei secular -omit this for homo
				bbb=1;

				for b=1:Dim
				for bb=1:Dim	

				delwbb=w_omega_D(1,b)-w_omega_D(1,bb); % Hetero nuclei secular -omit this for homo
				% calculate the relaxa element
				J2=J(bb,aa);
				J3=J(a,b);

				g1=0;
				k1=1;	


					% printf("J2:%0.10f,%0.10f   J3:%0.10f,%0.10f\n",real(J2),imag(J2),real(J3),imag(J3));


					Rel=0;
					for m1=-rank:rank
					Rel=Rel-J2*T1s(a,b,k1)*T2s(aa,bb,k1);

	     				%printf("RII:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(a,b,k1)),imag(T1s(a,b,k1)),real(T2s(aa,bb,k1)),imag(T2s(aa,bb,k1)),real(Rel),imag(Rel));



					Rel=Rel-J3*T1s(bb,aa,k1)*T2s(b,a,k1);

	    				%printf("RIII:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(bb,aa,k1)),imag(T1s(bb,aa,k1)),real(T2s(b,a,k1)),imag(T2s(b,a,k1)),real(Rel),imag(Rel));


					for g1=1:Dim

					if aa==bb
					Rel=Rel+J(g1,b)*T1s(a,g1,k1)*T2s(b,g1,k1);

				 %  printf("J1:%0.10f,%0.10f\n",real(J(g1,b)),imag(J(g1,b)));
	 			  % printf("RI:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(a,g1,k1)),imag(T1s(a,g1,k1)),real(T2s(b,g1,k1)),imag(T2s(b,g1,k1)),real(Rel),imag(Rel));



					end

					if a==b
					Rel=Rel+J(bb,g1)*T1s(g1,aa,k1)*T2s(g1,bb,k1);



				  % printf("J1:%0.10f,%0.10f\n",real(J(bb,g1)),imag(J(bb,g1)));
	 			  % printf("RIV:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(g1,aa,k1)),imag(T1s(g1,aa,k1)),real(T2s(g1,bb,k1)),imag(T2s(g1,bb,k1)),real(Rel),imag(Rel));


					end

					end %g1


					k1=k1+1;
					end %m1

				if Secular_approx ==1 
				if abs(delwaa-delwbb)<power(10,6) % Secular approx -omit for homo
				L_DC(aaa,bbb)=L_DC(aaa,bbb)+Rel;
				%printf("LDC[%d,%d,%d,%d]:%0.8f+%0.8f \n",a,aa,b,bb,real(L_DC(aaa,bbb)),imag(L_DC(aaa,bbb))); 
				end
				else
				L_DC(aaa,bbb)=L_DC(aaa,bbb)+Rel;
			%printf("LDC[%d,%d,%d,%d]:%0.8f+%0.8f \n",a,aa,b,bb,real(Rel),imag(Rel));  
				% printf("LR[%d,%d,%d,%d,%d,%d]:%0.8f+%0.8f \n",a,aa,b,bb,aaa,bbb,real(L_CC(aaa,bbb)),imag(L_CC(aaa,bbb)));  				pause
			%printf("LR[%d,%d,%d,%d]:%0.8f+%0.8f Cum:%0.8f+%0.8f \n",a,aa,b,bb,real(Rel),imag(Rel),real(L_CC(aaa,bbb)),imag(L_CC(aaa,bbb)));  		

				end % 	Secular_approx	


				bbb=bbb+1;
				end %bb
				end %b

				%printf("\n");				%pause
				aaa=aaa+1;
			%printf("\n");
			end %aa
			%printf("\n");

			end %a
			%printf("\n");
			%pause
			%L_DC
			% pause




		end %k

end %j
end %i


% Second cross interaction between CSA loop k vs dipolarloop ij
		for k=1:n_spins	% Dipolar

			xi2=real(xiC(k,k)); % CSA
			alphai=EA(k,1);
			betai=EA(k,2);
			gammai=EA(k,3);

			alphaij=alphai;
			betaij=betai;
			gammaij=gammai;

			K0 = 0.546274215;		% sqrt[15/(16*pi)]
			K1 = 1.092548431;		% sqrt[15/(4.0*pi)]
			K2 = 0.315391565;		% sqrt[5/(16.0*pi)]

			S2g = sin(2.0*gammaij);
			C2g = cos(2.0*gammaij);
			S2a = sin(2.0*alphaij);
			C2a = cos(2.0*alphaij);


			L = sin(betaij)*cos(alphaij);
			M = sin(betaij)*sin(alphaij);
			N = cos(betaij);


			etapart = (1.0+N*N)*C2a*C2g - 2.0*S2a*S2g*N;
			X = (L*L - M*M) + (eta/3.0)*etapart;
			Y = (3.0*N*N -1.0) + (eta*(1.0-N*N)*C2g);

			etapart = (eta/3.0)*((1+N*N)*S2a*S2g + 2.0*N*C2a*S2g);
			c1s(1,1) = K0 * (2.0*L*M + etapart);			% c -2 case

			etapart = (eta/3.0)*(M*S2g - L*N*C2g);
			c1s(1,2)= K1 * (L*N + etapart);				% c -1 case

			c1s(1,3)= K2 * (sqrt3*X*sin(chi) - Y*cos(chi));		% c  0 case

			etapart = (eta/3.0)*(L* + M*N*C2g);
			c1s(1,4)= K1 * (M*N - etapart);				% c +1 case

			c1s(1,5)= K2 * (sqrt3*X*cos(chi) + Y*sin(chi));		% c +2 case

			T1s=zeros(Dim,Dim,5);
			%{
			T1s(:,:,1)=H_V'*T_CS(:,:,1,k)*H_V; 
			T1s(:,:,2)=H_V'*T_CS(:,:,2,k)*H_V; 
			T1s(:,:,3)=H_V'*T_CS(:,:,3,k)*H_V; 
			T1s(:,:,4)=H_V'*T_CS(:,:,4,k)*H_V; 
			T1s(:,:,5)=H_V'*T_CS(:,:,5,k)*H_V; 
			%}

			T1s(:,:,1)=T_CS(:,:,1,k); 
			T1s(:,:,2)=T_CS(:,:,2,k); 
			T1s(:,:,3)=T_CS(:,:,3,k); 
			T1s(:,:,4)=T_CS(:,:,4,k); 
			T1s(:,:,5)=T_CS(:,:,5,k); 

			%{
			T1s(:,:,1)=T_CS(:,:,1,k); 
			T1s(:,:,2)=T_CS(:,:,2,k); 
			T1s(:,:,3)=T_CS(:,:,3,k); 
			T1s(:,:,4)=T_CS(:,:,4,k); 
			T1s(:,:,5)=T_CS(:,:,5,k); 
			%}

	for i=1:n_spins-1	% Dipolar
	for j=i+1:n_spins	% Dipolar


	xi1=xiD(i,j);	% Dipolar
	alphaij=Phis(i,j);
	betaij=Thetas(i,j);
	gammaij=0;

	K0 = 0.546274215;		% sqrt[15/(16*pi)]
	K1 = 1.092548431;		% sqrt[15/(4.0*pi)]
	K2 = 0.315391565;		% sqrt[5/(16.0*pi)]

	S2g = sin(2.0*gammaij);
	C2g = cos(2.0*gammaij);
	S2a = sin(2.0*alphaij);
	C2a = cos(2.0*alphaij);


	L = sin(betaij)*cos(alphaij);
	M = sin(betaij)*sin(alphaij);
	N = cos(betaij);


	etapart = (1.0+N*N)*C2a*C2g - 2.0*S2a*S2g*N;
	X = (L*L - M*M) + (eta/3.0)*etapart;
	Y = (3.0*N*N -1.0) + (eta*(1.0-N*N)*C2g);

	etapart = (eta/3.0)*((1+N*N)*S2a*S2g + 2.0*N*C2a*S2g);
	c2s(1,1) = K0 * (2.0*L*M + etapart);			% c -2 case

	etapart = (eta/3.0)*(M*S2g - L*N*C2g);
	c2s(1,2)= K1 * (L*N + etapart);				% c -1 case

	c2s(1,3)= K2 * (sqrt3*X*sin(chi) - Y*cos(chi));		% c  0 case

	etapart = (eta/3.0)*(L* + M*N*C2g);
	c2s(1,4)= K1 * (M*N - etapart);				% c +1 case

	c2s(1,5)= K2 * (sqrt3*X*cos(chi) + Y*sin(chi));		% c +2 case

	% c1s
	% pause


			T2s=zeros(Dim,Dim,5);

			%{
			T2s(:,:,1)=H_V'*T_D(:,:,1,i,j)*H_V; 
			T2s(:,:,2)=H_V'*T_D(:,:,2,i,j)*H_V; 
			T2s(:,:,3)=H_V'*T_D(:,:,3,i,j)*H_V; 
			T2s(:,:,4)=H_V'*T_D(:,:,4,i,j)*H_V; 
			T2s(:,:,5)=H_V'*T_D(:,:,5,i,j)*H_V; 
			%}

			T2s(:,:,1)=T_D(:,:,1,i,j); 
			T2s(:,:,2)=T_D(:,:,2,i,j); 
			T2s(:,:,3)=T_D(:,:,3,i,j); 
			T2s(:,:,4)=T_D(:,:,4,i,j); 
			T2s(:,:,5)=T_D(:,:,5,i,j); 

			%{
			T2s(:,:,1)=T_D(:,:,1,i,j); 
			T2s(:,:,2)=T_D(:,:,2,i,j); 
			T2s(:,:,3)=T_D(:,:,3,i,j); 
			T2s(:,:,4)=T_D(:,:,4,i,j); 
			T2s(:,:,5)=T_D(:,:,5,i,j); 
			%}


			
			% calculate the Spectral density matrix same for any mechanism depends only on CS and tau
			xi1xi2=xi1*xi2;
			for l1=1:Dim
			for l2=1:Dim	
				J(l1,l2)=0;

				wtr=real(w_omega_D(1,l1)-w_omega_D(1,l2));
	%printf("[%d,%d] omega[i]: %0.12f omega[j]: %0.12f w:  %0.12f\n",i,j,real(w_omega_D(1,l1)),real(w_omega_D(1,l2)),wtr);
 
				wtr=wtr*6.283185307;
				for m=1:5
				J_temp=Tau_eff(1,m)/(1+Tau_eff(1,m)*Tau_eff(1,m)*wtr*wtr);
				J(l1,l2)=J(l1,l2)+c1s(1,m)*c2s(1,m)*(J_temp);
	%printf("[%d] tau[i]: %0.12f omega: %0.12f c1[i]: %0.12f c2[i]: %0.12f J: %0.12f\n",m,Tau_eff(1,m),wtr,c1s(1,m),c1s(1,m),J(l1,l2));


				end
				J(l1,l2)=J(l1,l2)/5;
				%printf("J: %0.12f\n",J(l1,l2));
				%pause
				J(l1,l2)=J(l1,l2)*complex(xi1xi2,0);


			end %l1
			end %l2

			%{
			%printf("\n\n\nLoop[%d][%d][%d]:\n",k,i,j); pause
			%xi1xi2
			%T1s
			%T2s
			%c1s
			%%c2s
			%printf("\nJ\n");
			for l1=1:Dim
			for l2=1:Dim	
				%printf("%0.10f ",J(l1,l2));
			end
			%printf("\n");
			end
			%pause 			
			%printf("\n\n\n");
			%}


			%{
			printf("\n\n\nLoop[%d][%d][%d]:\n",i,j,k); pause
			xi1xi2
			T1s
			T2s
			c1s
			c2s
			printf("\nJ\n");
			for l1=1:Dim
			for l2=1:Dim	
				printf("%0.10f ",J(l1,l2));
			end
			printf("\n");
			end
			w_omega_D
			 pause 			
			%}
			% calculate the element of L_DC super operator
			aaa=1;


			Rel=0;
			J2=0;
			J3=0;


			for a=1:Dim
			for aa=1:Dim	
				delwaa=w_omega_D(1,a)-w_omega_D(1,aa); % Hetero nuclei secular -omit this for homo


				bbb=1;
				for b=1:Dim
				for bb=1:Dim	

				delwbb=w_omega_D(1,b)-w_omega_D(1,bb); % Hetero nuclei secular -omit this for homo
				% calculate the relaxa element
				J2=J(bb,aa);
				J3=J(a,b);

				g1=0;
				k1=1;	


					% printf("J2:%0.10f,%0.10f\n",real(J2),imag(J2));
					% printf("J3:%0.10f,%0.10f\n",real(J3),imag(J3));


					Rel=0;
					for m1=-rank:rank
					Rel=Rel-J2*T1s(a,b,k1)*T2s(aa,bb,k1);

	     				% printf("RII:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(a,b,k1)),imag(T1s(a,b,k1)),real(T2s(aa,bb,k1)),imag(T2s(aa,bb,k1)),real(Rel),imag(Rel));



					Rel=Rel-J3*T1s(bb,aa,k1)*T2s(b,a,k1);

	    				% printf("RIII:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(bb,aa,k1)),imag(T1s(bb,aa,k1)),real(T2s(b,a,k1)),imag(T2s(b,a,k1)),real(Rel),imag(Rel));


					for g1=1:Dim

					if aa==bb
					Rel=Rel+J(g1,b)*T1s(a,g1,k1)*T2s(b,g1,k1);

				  % printf("J1:%0.10f,%0.10f\n",real(J(g1,b)),imag(J(g1,b)));
	 			  % printf("RI:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(a,g1,k1)),imag(T1s(a,g1,k1)),real(T2s(b,g1,k1)),imag(T2s(b,g1,k1)),real(Rel),imag(Rel));



					end

					if a==b
					Rel=Rel+J(bb,g1)*T1s(g1,aa,k1)*T2s(g1,bb,k1);



				  % printf("J1:%0.10f,%0.10f\n",real(J(bb,g1)),imag(J(bb,g1)));
	 			 %  printf("RIV:T1s:%0.10f,%0.10f T2s:%0.10f,%0.10f Rel:%f,%f \n",real(T1s(g1,aa,k1)),imag(T1s(g1,aa,k1)),real(T2s(g1,bb,k1)),imag(T2s(g1,bb,k1)),real(Rel),imag(Rel));


					end

					end %g1


					k1=k1+1;
					end %m1

				if Secular_approx ==1 
				if abs(delwaa-delwbb)<power(10,6) % Secular approx -omit for homo
				L_DC(aaa,bbb)=L_DC(aaa,bbb)+Rel;
				%printf("LDC[%d,%d,%d,%d]:%0.8f+%0.8f \n",a,aa,b,bb,real(L_DC(aaa,bbb)),imag(L_DC(aaa,bbb))); 
				end
				else
				L_DC(aaa,bbb)=L_DC(aaa,bbb)+Rel;
				%printf("LDC[%d,%d,%d,%d]:%0.8f+%0.8f \n",a,aa,b,bb,real(Rel),imag(Rel));  
		

				end % 	Secular_approx
	
				bbb=bbb+1;
				end %bb
				end %b
				%printf("\n");				
				aaa=aaa+1;
			%printf("\n");
			end %aa
			%printf("\n");
			end %a
			%printf("\n");
			%pause
			%L_DC



		end %j
		end %i

end %k


% L_H0_EB=L_H0_EB-(L_DC+L_DD_RF+L_CC);	% Sign convention..Chemical exchange added (sign absorbed into kinetic matrix) but relaxation is substracted

L_Relax=L_DD_RF+L_CC-L_DC;

%}



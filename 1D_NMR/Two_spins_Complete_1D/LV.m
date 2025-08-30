% ==================================== 
% SIMULATION OF A Gradient SPECTRUM 
% Adopted from GAMMA
% ==================================== 
% NOTE: This program uses diprod.m, pauli.m and twospins.m. 
% All these files should be in the same folder. 


clear 
clc 
twospins


% INPUT PARAMETERS 
% ---------------- 

v1=25;
v2=50;
J12=10;
Dim=4;
MHZ=500;

Kex=500;	% Hz between spin 1-2


Gamma_1=267515255.000000;
Gamma_2=267515255.000000;	% 6.72669e+07; C13


RAD2HZ=1/(2*pi());
HZ2RAD=2*pi();
DEG2RAD=pi()/180;



% Setup the Hamiltonian


csham=-(v1*Izi + v2*Iiz);		% Chemical Shift Hamiltonian 
jham=J12*(Ixx+Iyy+Izz); 		% Scalar coupling Hamiltonian as in gamma 

H0=csham+jham; 			% Full Hamiltonian for evolution 
detop= Imi+Iim;			% Detection Operator


% Setup the Density matrix
sigma0=Izi+(Gamma_2/Gamma_1)*Iiz;
Sigma_L0=zeros(Dim*Dim,1);
% Flatten out the matrix to vector for initial denstiy matrix

%{
	m=1;
	for k=1:Dim
	for l=1:Dim
	Sigma_L(m,1)=sigma0(k,l);
	m=m+1;
	end
	end
%}

Sigma_L0=(reshape(transpose(sigma0),[],1));


% Apply 90 degree pulse

theta=pi/2; % Pulse angle 
Pulse_left = expm(-complex(0,1)*(Iyi+Iiy)*theta);


%{
Pulse_left = expm(-complex(0,1)*(Iyi+Iiy)*theta);
Pulse_right = expm(complex(0,1)*(Iyi+Iiy)*theta);
sigmap_y=Pulse_left*sigma0*Pulse_right;	% First pulse

% L_pulse_90_H12y=diprod(expm(-complex(0,1)*(Iyi+Iiy)*theta),conj(expm(-complex(0,1)*(Iyi+Iiy)*theta)));
%}

L_pulse_90_H12y=diprod(Pulse_left,conj(Pulse_left));
sigmap_y=L_pulse_90_H12y*Sigma_L0;


Sigma_L=zeros(Dim*Dim,1);
Sigma_L=sigmap_y;




% Liouville parameters
DIM=Dim*Dim;
L_i=eye(Dim,Dim);
L_I=eye(DIM,DIM);
L_H0=-(diprod(H0,L_i)-diprod(L_i,H0));



% Vectorize Sigma and detector as superket

Detect_L=zeros(1,Dim*Dim);


% Flatten out the matrix to vector for detector
	%{
	m=1;
	for k=1:Dim
	for l=1:Dim
	Detect_L(1,m)=detop(k,l);
	m=m+1;
	end
	end
	%}
	
Detect_L=transpose(reshape(transpose(detop),[],1));



% FID parameters
%---------------
NyqF=100;
SW=2*NyqF;
npts = 2048;
td = 1/SW;
freq=-SW/2:SW/(npts-1):SW/2;


T2=2;
T_relax = exp(-td/T2);



FID=zeros(1,npts);
Time=zeros(1,npts);
U_left=eye(DIM,DIM);
for i=1:npts

	Time(1,i)=(i-1)*td;
	if i==1
	sigmap_L=Sigma_L;
	end


	sigmap_L=(U_left*sigmap_L).*T_relax; %.*T_relax;	%.*T_relax;
	FID(1,i)=(conj(Detect_L)*sigmap_L);		% Detection 
	U_left = expm(-complex(0,1)*2*pi*L_H0*td);



end



% plot
%-----

Ft=fftshift(fft(FID)); % Fourier Transform 
figure();	
plot(Time(1,:),FID(1,:),'b-','LineWidth',1.5);
figure();	
plot(freq(1,:),Ft(1,:),'b-','LineWidth',1.5);

outfile=fopen("Results_LV","w");
for i=1:npts
fprintf(outfile,"FID: %f  %f  %f FT: %f  %f  %f\n",Time(1,i),real(FID(1,i)),imag(FID(1,i)),freq(1,i),real(Ft(1,i)),imag(Ft(1,i)));
end
fclose(outfile);


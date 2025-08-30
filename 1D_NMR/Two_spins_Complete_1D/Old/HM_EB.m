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



Gamma_1=267515255.000000;
Gamma_2=267515255.000000;	% 6.72669e+07; C13


RAD2HZ=1/(2*pi());
HZ2RAD=2*pi();
DEG2RAD=pi()/180;



% Setup the Hamiltonian

sigma0=Izi+(Gamma_2/Gamma_1)*Iiz;

csham=-(v1*Izi + v2*Iiz);		% Chemical Shift Hamiltonian 
jham=J12*(Ixx+Iyy+Izz); 		% Scalar coupling Hamiltonian as in gamma 

H0=csham+jham; 			% Full Hamiltonian for evolution 

[H_V,H_D]=eig(H0);


detop= Imi+Iim;			% Detection Operator

detop_EB=H_V'*detop*H_V;

% Apply 90 degree pulse

theta=pi/2; % Pulse angle 

Pulse_left = expm(-complex(0,1)*(Iyi+Iiy)*theta);
Pulse_right = expm(complex(0,1)*(Iyi+Iiy)*theta);

%sigmap_y=Pulse_left*sigma0*Pulse_right;	% First pulse
sigma0_EB=H_V'*sigma0*H_V;
sigmap_y_EB=(H_V'*Pulse_left*H_V*sigma0_EB*H_V'*Pulse_right*H_V);
%sigmap_y_EB=(H_V'*Pulse_left*H_V*sigma0*H_V'*Pulse_right*H_V);
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
U_left=eye(Dim,Dim);
U_right=eye(Dim,Dim);
for i=1:npts

	Time(1,i)=(i-1)*td;
	if i==1
	sigmap=sigmap_y_EB;
	end


	sigmap=(U_left*sigmap*U_right).*T_relax; %.*T_relax;	%.*T_relax;
	FID(1,i)=trace(detop_EB*sigmap);		% Detection 
	U_left = expm(-1*complex(0,1)*2*pi*H_D*td);
	U_right = expm(1*complex(0,1)*2*pi*H_D*td);
	detop_EB*sigmap
	FID(1,i)
	pause

end



% plot
%-----

Ft=fftshift(fft(FID)); % Fourier Transform 
%Ft=(fft(FID));
figure();	
plot(Time(1,:),FID(1,:),'b-','LineWidth',1.5);
figure();	
plot(freq(1,:),Ft(1,:),'b-','LineWidth',1.5);

outfile=fopen("Results_HM_EB","w");
for i=1:npts
fprintf(outfile,"FID: %f  %f  %f FT: %f  %f  %f\n",Time(1,i),real(FID(1,i)),imag(FID(1,i)),freq(1,i),real(Ft(1,i)),imag(Ft(1,i)));
end
fclose(outfile);


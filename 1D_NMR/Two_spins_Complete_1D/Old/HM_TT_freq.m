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
MHZ=500;
n_spins=2;
Dim=power(n_spins,2);	% depends on the spin states of the nuclei alpha,beta =2 for spin 1/2.


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

sigmap_y_EB=(H_V'*Pulse_left*H_V*sigma0*H_V'*Pulse_right*H_V);
% FID parameters
%---------------
NyqF=100;
SW=2*NyqF;
npts = 2048;
td = 1/SW;
freq=-SW/2:SW/(npts-1):SW/2;


T2=2;
T_relax = exp(-td/T2);


Ft=zeros(1,npts);
Time=zeros(1,npts);
U_left=eye(Dim,Dim);
U_right=eye(Dim,Dim);

U_left = expm(-1*complex(0,1)*2*pi*H_D*td);
U_right = expm(1*complex(0,1)*2*pi*H_D*td);



% Generate the full transition table using H eigen basis element rather than Propogator (see HM_TT_time_alt.m)
Ttable=zeros(Dim*Dim,3); 	% A=Intensity (I), B= transition frequency (w), C= Relaxation (R)

T_index=1;
for i=1:Dim
for j=1:Dim
	Ttable(T_index,1)=detop_EB(i,j);
	Ttable(T_index,2)=-H_D(j,j)+H_D(i,i);
	Ttable(T_index,3)=sigmap_y_EB(j,i);
	T_index=T_index+1;
end % col index
end % row index


% Remove enteries whose enteries A (arising from Detector matrix), C (from sigma0) are zeros
Pos=1;
T_index=1;
Ttable_final=zeros(Dim*Dim,3);

for i=1:Dim
for j=1:Dim

	if (Ttable(T_index,1) && Ttable(T_index,3)) ~= 0
	Ttable_final(Pos,1)=detop_EB(i,j);
	Ttable_final(Pos,2)=-H_D(j,j)+H_D(i,i);
	Ttable_final(Pos,3)=sigmap_y_EB(j,i);
	Pos=Pos+1;
	end
	T_index=T_index+1;
end % col index
end % row index


Time=(0:td:td*(npts-1));
% FID using the reduced Ttable : Use this FID loop for faster calculation with less transition enteries
%{
for i=1:Pos
for j=1:npts
	FID(1,j)=FID(1,j)+(Ttable_final(i,1)*power(Ttable_final(i,2),j-1)*Ttable_final(i,3))*power(T_relax,j-1);
end % j (fid)
end % i (transition)
%}

% FID using the full Ttable just to cross check the exactness with the other direct methods

%LW=0.5*2*pi;%Hz->w
T2=2;	% sec
Norm_factor=npts/sqrt(2*pi);
for i=1:Dim*Dim
j=1;
for freq_1=-SW/2:SW/(npts-1):SW/2

	%B=complex(0,1)*2*pi*(freq_1-Ttable(i,2));

	%Ft(1,j)=Ft(1,j)+Ttable(i,1)*(Ttable(i,3)-2*pi*complex(0,1)*(freq_1-Ttable(i,2)))/(power(Ttable(i,3),2)+power(2*pi*(freq_1-Ttable(i,2)),2)); % close

	B=complex(0,1)*2*pi*(freq_1-Ttable(i,2))-2*pi/T2;
	Ft(1,j)=Ft(1,j)+Norm_factor*Ttable(i,1)*Ttable(i,3)*(-real(B)+complex(0,1)*imag(B))/(power(real(B),2)+power(imag(B),2)); 

%{
	B=2*pi*(freq_1-Ttable(i,2));
	Ft(1,j)=Ft(1,j)+Ttable(i,1)*Ttable(i,3)*(-imag(B)-complex(0,1)*(real(B)+LW))/(power(real(B)+LW,2)+power(imag(B),2)); 

%}

	j=j+1;
end % j (fid)
end % i (transition)



% plot
%-----

%Ft=fftshift(fft(FID)); % Fourier Transform 
%Ft=(fft(FID));
%figure();	
%plot(Time(1,:),FID(1,:),'b-','LineWidth',1.5);
figure();	
plot(freq(1,:),(Ft(1,:)),'b-','LineWidth',1.5);
figure();
plot(Ft(1,:),'b-','LineWidth',1.5);
outfile=fopen("Results_HM_TT_freq","w");
for i=1:npts
fprintf(outfile,"FT: %f  %f  %f\n",freq(1,i),real(Ft(1,i)),imag(Ft(1,i)));
end
fclose(outfile);


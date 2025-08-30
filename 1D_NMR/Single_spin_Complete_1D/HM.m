% ==================================== 
% SIMULATION OF A Gradient SPECTRUM 
% Adopted from GAMMA
% ==================================== 
% NOTE: This program uses diprod.m, pauli.m and twospins.m. 
% All these files should be in the same folder. 


clear 
clc 
nspins

mkdir('./Results');

% INPUT PARAMETERS 
% ---------------- 

v1=25;
Dim=2;
MHZ=500;



Gamma_1=267515255.000000;

RAD2HZ=1/(2*pi());
HZ2RAD=2*pi();
DEG2RAD=pi()/180;



% Setup the Hamiltonian

sigma0=Iz;

csham=-(v1*Iz);		% Chemical Shift Hamiltonian 

H0=csham; 			% Full Hamiltonian for evolution 
detop= Im;			% Detection Operator


% Apply 90 degree pulse

theta=pi/2; % Pulse angle 

Pulse_left = expm(-complex(0,1)*(Iy)*theta);
Pulse_right = expm(complex(0,1)*(Iy)*theta);

sigmap_y=Pulse_left*sigma0*Pulse_right;	% First pulse


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
	sigmap=sigmap_y;
	end


	sigmap=(U_left*sigmap*U_right).*T_relax; %.*T_relax;	%.*T_relax;
	FID(1,i)=trace(detop*sigmap);		% Detection 
	U_left = expm(-1*complex(0,1)*2*pi*H0*td);
	U_right = inv(U_left);

end



% plot
%-----

Ft=fftshift(fft(FID)); % Fourier Transform 
%Ft=(fft(FID));
figure();	
plot(Time(1,:),FID(1,:),'b-','LineWidth',1.5);
saveas(gcf, sprintf('1D(FID)'), 'pdf');
figure();	
plot(freq(1,:),Ft(1,:),'b-','LineWidth',1.5);
saveas(gcf, sprintf('1D(FT)'), 'pdf');


outfile=fopen("./Results/Results_HM","w");
for i=1:npts
fprintf(outfile,"FID: %f  %f  %f FT: %f  %f  %f\n",Time(1,i),real(FID(1,i)),imag(FID(1,i)),freq(1,i),real(Ft(1,i)),imag(Ft(1,i)));
end
fclose(outfile);




outfile=fopen("./Results/Latex_fid_data","w");
fprintf(outfile,"[{");
for i=1:npts
if i<npts
fprintf(outfile,"{%f, %f},\n",Time(1,i),real(FID(1,i)));
else
fprintf(outfile,"{%f, %f}}]",Time(1,i),real(FID(1,i)));
end
end
fclose(outfile);

outfile=fopen("./Results/Latex_ft_data","w");
fprintf(outfile,"[{");
for i=1:npts
if i<npts
fprintf(outfile,"{%f, %f},\n",freq(1,i),real(Ft(1,i)));
else
fprintf(outfile,"{%f, %f}}]",freq(1,i),real(Ft(1,i)));
end
end
fclose(outfile);
%}

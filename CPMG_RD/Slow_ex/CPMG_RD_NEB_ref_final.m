
% SIMULATION OF A TWO-SPIN 2D HSQC SPECTRUM with Relaxation+Chemical exchange containing CPMG Relaxation dispersion element


% Std HSQC:

%	I: 90x--tau--180x--tau--90y--|  	|-------90x-------tau--180x--tau--t2
%			         	  CPMG
%	S: -----tau--180x--tau--90x--|  	|---t1--90q[x,y]--tau--180x--tau

% CPMG  : Tollinger

% ---------------|--Tau_RD--------Tau_RD--|   ----90y180x90y----|--Tau_RD--------Tau_RD--|-------------tau--180x--tau--
% ---------------|--Tau_RD--180y--Tau_RD--|[n]-------180x-------|--Tau_RD--180x--Tau_RD--|[n]----------tau--180x--tau--

%<--------------------IzSy----------------------------><--------------------Sx----------------------><----- IzSy------>

% Ref is wiithout the cpmg element


% Fast exchange:
% --------------
 
% Change pop (95:5), del W (0.1), Kex 10000, P element :I 90x180y90x between CPMG blocks TR=40ms
% [N_factor=CPMG_N=4; TR=44.44ms; N_pulses=2,4,8,16,32,64 [or in Hz :45 Hz, 90, 180, 360, 720, 1440 Hz]






% by Jana.K
% Ackowledgement:  Scott A Smith, T. S. Mahesh, Ilya Kuprov's, Tata Gopinath
% ===========================================================================

clear 
clc 
tic 
nspins
npulses

% CPMG RD specific inputs
TR_SCALE_PERIOD=4; 	% an integer say 2 or 4 or 5  etc, which extends the TR period
N_CPMG_PULSES=[0]; 	%  Each entry correspond to separate ext with the entry number as the 
			%  number of pi pulses applied during TR/2 period. its a row vector,for 3 expts: e.g [2 4 8]
			% for 10 expts =[2 4 6 8 12 16 20 24 32 64]; 
N_expts=size(N_CPMG_PULSES); % take the col value for looping for different expts, N_expts(1,2)

% INPUT PARAMETERS 
% ---------------- 

v_iso_ppm_orig_1=[6.0 102.0];	% Proton, Nitrogen isotropic CS (ppm) of first NH spin system 
v_iso_ppm_orig_2=[7.0 103.0];	% Proton, Nitrogen isotropic CS (ppm) of second NH spin system 

offset_ppm=[6.0 102.0];

v_iso_ppm_1=v_iso_ppm_orig_1-offset_ppm; 	% Proton, Nitrogen isotropic CS (ppm) of first NH spin system 
						% with respect to offset/center of Spectrum.
v_iso_ppm_2=v_iso_ppm_orig_2-offset_ppm;


% Population
%-----------
Population_spin_Eq=[1 1]; % The proportion of H, N nuclei

% Proportion of exchange species A<===>B of HN
Pop_A=[95]; % The number of enteries here detemine the number of titrations 
Pop_B=100-Pop_A;
K=power(10,1); % set value as Fast: 100000, Intermediate: 500, Slow: 10 % Octave expm() hv issues with K>10K matlab doesnt
		% used matlab here.. Kex or conformational ex, (k_on+k_off)
		
Titration_count=size(Pop_A);
Titration_no=Titration_count(1,2);

Population_Ex=[Pop_A Pop_B];
min_Pop=min(Population_Ex);
Population_norm=Population_Ex./min_Pop;



% Spin quantum number
%-------------------
S_1=[1/2 1/2]; % Spin quantum number 1H,15N
S_2=[1/2 1/2]; % Spin quantum number 1H,15N

J12=90;
J34=90;
J_I=[0 J12 ; J12 0];
J_II=[0 J34; J34 0]; % J couple matrix


J_HN=mean([J12 J34]);
Tau=1/(4*J_HN); 


n_spins_1=2;
n_spins_2=2;

Dim=[power(n_spins_1,2) power(n_spins_2,2)];

MHZ_H=500;	% Spectrometer field strength

Dim_1=S_1.*2+1;
Dim_2=S_2.*2+1;

Dim_T_1=prod(Dim_1(1,:));
Dim_T_2=prod(Dim_2(1,:));

Dim_T=Dim_T_1+Dim_T_2;


% Calculating the MHZ for H,15N for the desired MHZ with respecto values at 400 MHz

Rfreq_1=[400.130 -40.561];	% Ref spectrometer 400 MHz H,15N
Rfreq_2=[400.130 -40.561];

Gamma_H=2.67515255*power(10,8);
Rfreq_H=400.130;

Gamma_1=Rfreq_1.*Gamma_H/Rfreq_H; % Gamma calculated based on Hetero nuclei freq at 400 MHz and gamma of H
MHZ_1=Gamma_1.*(MHZ_H/Gamma_H);
v_Hz_1=(diag(v_iso_ppm_1*eye(n_spins_1,n_spins_1))*MHZ_1')';

Gamma_2=Rfreq_2.*Gamma_H/Rfreq_H; % Gamma calculated based on Hetero nuclei freq at 400 MHz and gamma of H
MHZ_2=Gamma_2.*(MHZ_H/Gamma_H);
v_Hz_2=(diag(v_iso_ppm_2*eye(n_spins_2,n_spins_2))*MHZ_2')';


swh1=4000; % Spectral width for F1 Dimension 
swh2=4000; % Spectral width for F2 Dimension 
td1=256; % t1 Time domain size (Indirect dimension)
td2=1024; % t2 Time domain size  (Direct dimension)

% CALCULATING TIME AND FREQUENCY AXIS 
% ----------------------------------- 

dw1=1/swh1; % Dwell time 
t1=0:dw1:(td1-1)*dw1; % Time axis 
f1=-swh1/2:swh1/(td1-1):swh1/2;

f1=f1./MHZ_1(1,2);

dw2=1/swh2; % Dwell time 
t2=0:dw2:(td2-1)*dw2; % Time axis 
f2=-swh2/2:swh2/(td2-1):swh2/2;
f2=f2./MHZ_1(1,1);

% INITIAL DENSITY MATRIX AND HAMILTONIANS 
% --------------------------------------- 


H1_csham=-(v_Hz_1(1,1)*HzNi); % Chemical Shift Hamiltonian 
N1_csham=-(v_Hz_1(1,2)*HiNz); % Chemical Shift Hamiltonian 
H2_csham=-(v_Hz_2(1,1)*HzNi); % Chemical Shift Hamiltonian 
N2_csham=-(v_Hz_2(1,2)*HiNz); % Chemical Shift Hamiltonian 

T_csham_1=H1_csham+N1_csham; 	% Chemical Shift Hamiltonian of first component or spin sys's NH (I) 
T_jham_1=J_I(1,2)*HzNz;		% Scalar coupling Hamiltonian of first component or spin sys's NH (I) 
H0_1=T_csham_1+T_jham_1;	% Total Hamiltonian of first component or spin sys's NH (I) 
H0_1_decp=T_csham_1;



T_csham_2=H2_csham+N2_csham; 	% Chemical Shift Hamiltonian of first component or spin sys's NH (II) 
T_jham_2=J_II(1,2)*HzNz;	% Scalar coupling Hamiltonian of first component or spin sys's NH (II) 
H0_2=T_csham_2+T_jham_2;	% Total Hamiltonian of first component or spin sys's NH (II) 
H0_2_decp=T_csham_2;







[H_V1,H_D1]=eig(H0_1);
H0_1_EB=H_D1;

[H_V1_decp,H_D1_decp]=eig(H0_1_decp);
H0_1_decp_EB=H_D1_decp;

[H_V2,H_D2]=eig(H0_2);
H0_2_EB=H_D2;

[H_V2_decp,H_D2_decp]=eig(H0_2_decp);
H0_2_decp_EB=H_D2_decp;



% Initial density matrix
%-----------------------
Comp_I=[1 0; 0 0];
Comp_II=[0 0; 0 1];

Sigma_Comp_I=[1; 0];
Sigma_Comp_II=[0; 1];



sigma0_1=Population_Ex(1,1)*Population_spin_Eq(1,1)*((Gamma_1(1,1)/Gamma_1(1,1))*HzNi+(Gamma_1(1,2)/Gamma_1(1,1))*HiNz);
sigma0_2=Population_Ex(1,2)*Population_spin_Eq(1,2)*((Gamma_2(1,1)/Gamma_2(1,1))*HzNi+(Gamma_2(1,2)/Gamma_2(1,1))*HiNz);




sigma0_Comp_1=diprod(Sigma_Comp_I,sigma0_1);
sigma0_Comp_2=diprod(Sigma_Comp_II,sigma0_2);

sigma0=sigma0_Comp_1+sigma0_Comp_2;

sigma0_L=(reshape(transpose(sigma0),[],1));
sigma0_L1=(reshape(transpose(sigma0_Comp_1),[],1));
sigma0_L2=(reshape(transpose(sigma0_Comp_2),[],1));

%%rho_L_1=transpose(reshape(transpose(rho_1),[],1));	
%rho_L_1=(U_left_tau*rho_L_1); %.*T_relax;
%rho_1=reshape(C,size(sigma0)')' 



% Detector
%---------

detop_Comp_1=diprod(Sigma_Comp_I,HmNi);
detop_Comp_2=diprod(Sigma_Comp_II,HmNi);


detop=detop_Comp_1+detop_Comp_2;

detop_L=transpose(reshape(transpose(detop),[],1));


% Based on Ilya's suggestions generate propogators outside for loop
%----------------------------
fid_cos_1=zeros(td1,td2);
fid_cos_2=zeros(td1,td2);
fid_sin_1=zeros(td1,td2);
fid_sin_2=zeros(td1,td2);

fid_cos=zeros(td1,td2);
fid_sin=zeros(td1,td2);

ft_1=zeros(td1,td2);
ft_2=zeros(td1,td2);
ft_3=zeros(td1,td2);
ft_4=zeros(td1,td2);
fid=zeros(td1,td2);
ft=zeros(td1,td2);




% Pulses along X,Y for H and N and its propogators
%--------------------------------------------------

% H
%---

phi=pi;		% 0: along 0: a, pi/2:-b, pi:-a, -pi/2:+b, -pi=0:a ; for arbitrary angle starting from y axis pulse
		% if a=x then b=y and viceversa

% Pulses in Liouville space




% Test print the pulses
%----------------------


% Precalculate Propogator for t2 evolution
% Note the prop for t1 delay alone needs sign alteration to get correct CS in the indirect dimension, which seems to be correct

% Liouville parameters
%---------------------
DIM=Dim_T_1*Dim_T_1+Dim_T_2*Dim_T_2;


L_i=eye(Dim_T_1,Dim_T_1);
L_H0_1=complex(0,1)*2*pi*(diprod(H0_1,L_i)-diprod(L_i,H0_1));
L_H0_1=diprod(Comp_I,L_H0_1);

L_H0_1_decp=-complex(0,1)*2*pi*(diprod(H0_1_decp,L_i)-diprod(L_i,H0_1_decp)); % for proper CS -ve needed
L_H0_1_decp=diprod(Comp_I,L_H0_1_decp);

L_i=eye(Dim_T_2,Dim_T_2);
L_H0_2=complex(0,1)*2*pi*(diprod(H0_2,L_i)-diprod(L_i,H0_2));
L_H0_2=diprod(Comp_II,L_H0_2);

L_H0_2_decp=-complex(0,1)*2*pi*(diprod(H0_2_decp,L_i)-diprod(L_i,H0_2_decp)); % for proper CS -ve needed
L_H0_2_decp=diprod(Comp_II,L_H0_2_decp);

L_H0=L_H0_1+L_H0_2;
L_H0_decp=L_H0_1_decp+L_H0_2_decp;

% all H decouple and N : not used
L_i=eye(Dim_T_2,Dim_T_2);
L_H0_H1=complex(0,1)*2*pi*(diprod(H1_csham,L_i)-diprod(L_i,H1_csham));
L_H0_H1=diprod(Comp_I,L_H0_H1);

L_H0_H2=complex(0,1)*2*pi*(diprod(H2_csham,L_i)-diprod(L_i,H2_csham));
L_H0_H2=diprod(Comp_II,L_H0_H2);


L_H0_N1=complex(0,1)*2*pi*(diprod(N1_csham,L_i)-diprod(L_i,N1_csham));
L_H0_N1=diprod(Comp_I,L_H0_N1);

L_H0_N2=complex(0,1)*2*pi*(diprod(N2_csham,L_i)-diprod(L_i,N2_csham));
L_H0_N2=diprod(Comp_II,L_H0_N2);






% Calculate Relaxation matrix 

	% Common properties
	Secular_approx=1;	% Relaxation with secular approximation (only for hetero nuclei: 1; homo set: 0)
	% 0,0 gives correct results in matlab but others freq are correctly represented by nmrpipe	
	Tau_D=[100 100 100].*power(10,-9); 	% correlation along x,y,z in nano sec


	
	% Specific for 1st H-N pair
	n_spins=n_spins_1;	
	J_H1N2 = J12; % H-N coupling in Hz
	S=S_1; % Spin quantum number
	v_ppm=v_iso_ppm_1; 	% Proton ppm	 Nitrogen ppm

	% CSA Tensor
	%---------- 
	v_iso=v_iso_ppm_1; %[0.0 0.0];	% Isotropic PPM
	v_delzz=[-16 -160];
	v_eta=[0 0];

	Coord_1=[0.9, -0.3, -0.4];	% center of spin 1
	Coord_2=[-0.2, 1.0, -0.4];	% center of spin 2

	EA_1=[90.0, 0.0, 30.0];		% Alpha, Beta, Gamma : orientation of spin 1
	EA_2=[0.0, 90.0, 20.0];		% Alpha, Beta, Gamma : orientation of spin 2

%REX_1=Relax_matrix(v_ppm,Secular_approx,J_H1N2,n_spins,MHZ_H,v_iso,v_delzz,v_eta,Coord_1,Coord_2,EA_1,EA_2,S,Tau_D);
[L_DD,L_CC,L_DC]=Relax_matrix(v_ppm,Secular_approx,J_H1N2,n_spins,MHZ_H,v_iso,v_delzz,v_eta,Coord_1,Coord_2,EA_1,EA_2,S,Tau_D);

%REX_1=L_DD+L_CC+L_DC; %-ve sign for L_DC will mirror the relaxed couplets from left to right
REX_1=L_DD;	

	% Specific for 1st H-N pair
	n_spins=n_spins_2;	
	J_H1N2 = J34; % H-N coupling in Hz
	S=S_2; % Spin quantum number
	v_ppm=v_iso_ppm_2; 	% Proton ppm	 Nitrogen ppm

	% CSA Tensor
	%---------- 
	v_iso=v_iso_ppm_2;% [1.0 1.0];	% Isotropic PPM
	v_delzz=[-16 -160];
	v_eta=[0 0];

	Coord_1=[3.9, 2.7, 2.6];	% center of spin 1  +3 translation
	Coord_2=[2.8, 4.0, 2.6];	% center of spin 2  +3 translation

	EA_1=[90.0, 0.0, 30.0];		% Alpha, Beta, Gamma : orientation of spin 1
	EA_2=[0.0, 90.0, 20.0];		% Alpha, Beta, Gamma : orientation of spin 2

[L_DD,L_CC,L_DC]=Relax_matrix(v_ppm,Secular_approx,J_H1N2,n_spins,MHZ_H,v_iso,v_delzz,v_eta,Coord_1,Coord_2,EA_1,EA_2,S,Tau_D);

%REX_2=L_DD+L_CC+L_DC; %-ve sign for L_DC will mirror the relaxed couplets from left to right
REX_2=L_DD;


L_i2=eye(Dim_T_2/2,Dim_T_2/2);

REX_Comp_1=diprod(Comp_I,REX_1);
REX_Comp_2=diprod(Comp_II,REX_2);
REX=REX_Comp_1+REX_Comp_2;		% Decoupling along both H+N axis (only H,N CS evolves)



% Chemical exchange and relaxation parameters

p_a=Population_norm(1,1);
p_b=Population_norm(1,2);
hs1=1;
hs2=1;
hs3=1;
hs4=1;

POP=[p_a p_b];
hs= [hs1 hs2 hs3 hs4];

KEX=Kex_matrix(K,POP,hs);

L_H0_decp_final=L_H0_decp+REX;
L_H0_cpmg=L_H0+REX+KEX;
L_H0_nocpmg=L_H0+REX;

L_H0_decp_H=L_H0_H1+L_H0_H2+REX+KEX;
L_H0_decp_N=L_H0_N1+L_H0_N2+REX+KEX;


% Propagator for delays

U_left_tau = expm(-L_H0_nocpmg*Tau); 	% expm(-L_H0_final*Tau); other than cpmg, t1, t2
U_left_t1 = expm(-L_H0_decp_final*(dw1/2)); % for correct CS sign is reversed for H0_decp_final i.e expm(+)
U_left_t2 = expm(-L_H0_decp_final*dw2);	% L_H0_decp_H ( bcs in expt 180 pulse in middle of t1 will reverse the sign of coherence)

% CPMG RD parameters	
CPMG_N=TR_SCALE_PERIOD; % an integer say 4,5 etc
TR=(CPMG_N/J_HN);	% Total relaxation period

	% NMR Pipe script
	mkdir('./NMRPipe');
	File_Result=sprintf('./NMRPipe/NMRpipe_script_ref');
	outfile_2 = fopen(File_Result,"w");	
	fprintf(outfile_2,'#!/bin/tcsh\n');	
	fprintf(outfile_2,'\n# make this file executable by typing in a cmd window or terminal: chmod 777 ./NMRpipe_script');	
	fprintf(outfile_2,'\n# adjust the phases P0,P1 in both dim H,N, if necessary \n');	

for cpmg_expt=1:N_expts(1,2)

	% CPMG RD parameters
	N_CPMG_Pulse=N_CPMG_PULSES(1,cpmg_expt); 	% pulses in half block (TR/2) ,it determines 
	Tau_TR=(TR/2)/(2*N_CPMG_Pulse); 		% (one half) of the spin echo period for each pi pulse applied
							% the field strength during TR by:2, 4 ,16
							
							
	%[L_H0_cpmg_v,L_H0_cpmg_d]=eig(L_H0_cpmg);						
	%U_left_cpmg_diag_tr_tau = expm(-L_H0_cpmg_d*Tau_TR);						
	U_left_cpmg_tr_tau = expm(-L_H0_cpmg*Tau_TR);	
	U_left_cpmg_tau = expm(-L_H0_cpmg*Tau);		
	
	
	

	File_Result=sprintf('./NMRPipe/CPMG_RD_%dPulse_%.0f_kex_%.0fA-%.0fB_Pop.asc',N_CPMG_Pulse,K,Pop_A,Pop_B);
	outfile_1 = fopen(File_Result,"w");



for n=1:td1
% Cos modulation					% t1

    	rho_L1= sigma0_L; 
     	rho_L2= sigma0_L; 



   	rho_L1=(L_pulse_90_N1y)*rho_L1;	
   	rho_L2=(L_pulse_90_N2y)*rho_L2;	

   	     	
	% 90 H pulse +x
        rho_L1= L_pulse_90_H1x*rho_L1;			
        rho_L2= L_pulse_90_H2x*rho_L2;		

	% Tau evolution -1
	rho_L1=(U_left_tau*rho_L1); %.*T_relax;	
	rho_L2=(U_left_tau*rho_L2); %.*T_relax;	
   	
 	% 180 H pulse +x and 180 N pulse +x 
   	rho_L1=(L_pulse_180_H1x * L_pulse_180_N1x)*rho_L1;		
   	rho_L2=(L_pulse_180_H2x * L_pulse_180_N2x)*rho_L2;

	% Tau evolution -2
	rho_L1=(U_left_tau*rho_L1); %.*T_relax;	
	rho_L2=(U_left_tau*rho_L2); %.*T_relax;
		
 	% 90 H pulse +y and 90 N pulse +x 
   	rho_L1=(L_pulse_90_H1y)*rho_L1;	
   	rho_L2=(L_pulse_90_H2y)*rho_L2;
   	
   	rho_L1=(L_pulse_90_N1x)*rho_L1;	
   	rho_L2=(L_pulse_90_N2x)*rho_L2;

%{		
		% CPMG RD element left half
		for cpmg_loop=1:N_CPMG_Pulse
			rho_L1=U_left_cpmg_tr_tau*rho_L1; 	% T_relax - tau;	
			rho_L2=U_left_cpmg_tr_tau*rho_L2; 	% T_relax - tau;			
			
		   	rho_L1=L_pulse_180_N1y*rho_L1; 	% pi pulse y
		   	rho_L2=L_pulse_180_N2y*rho_L2; 	% pi pulse y
		   	
		   	rho_L1=U_left_cpmg_tr_tau*rho_L1; 	% T_relax- tau;	
		   	rho_L2=U_left_cpmg_tr_tau*rho_L2; 	% T_relax- tau;
		   	
	   		
		end

		%  CPMG : P element 
			rho_L1=U_left_cpmg_tau*rho_L1; 	% T_relax - tau;	
			rho_L2=U_left_cpmg_tau*rho_L2; 	% T_relax - tau;

		   	rho_L1=(L_pulse_90_H1y*L_pulse_180_H1x*L_pulse_90_H1y)*rho_L1; 	% p element pulses H		
		   	rho_L2=(L_pulse_90_H2y*L_pulse_180_H2x*L_pulse_90_H2y)*rho_L2; 	% p element pulses H	

		   	rho_L1=(L_pulse_180_N1x)*rho_L1; 	% p element pulses N			   		
		   	rho_L2=(L_pulse_180_N2x)*rho_L2; 	% p element pulses N	

			rho_L1=U_left_cpmg_tau*rho_L1; 	% T_relax - tau;		   		
			rho_L2=U_left_cpmg_tau*rho_L2; 	% T_relax - tau;	

			
		% CPMG RD element right half
		for cpmg_loop=1:N_CPMG_Pulse
			rho_L1=U_left_cpmg_tr_tau*rho_L1; 	% T_relax - tau;	
			rho_L2=U_left_cpmg_tr_tau*rho_L2; 	% T_relax - tau;	

		   	rho_L1=L_pulse_180_N1x*rho_L1; 	% pi pulse x
		   	rho_L2=L_pulse_180_N2x*rho_L2; 	% pi pulse x

		   	rho_L1=U_left_cpmg_tr_tau*rho_L1; 	% T_relax- tau;	
		   	rho_L2=U_left_cpmg_tr_tau*rho_L2; 	% T_relax- tau;	
		end


	% To recover IzSy from Sx

	% Tau evolution -1
	rho_L1=(U_left_tau*rho_L1); %.*T_relax;	
	rho_L2=(U_left_tau*rho_L2); %.*T_relax;	
   	
 	% 180 H pulse +x and 180 N pulse +x 
   	rho_L1=(L_pulse_180_H1x * L_pulse_180_N1x)*rho_L1;		
   	rho_L2=(L_pulse_180_H2x * L_pulse_180_N2x)*rho_L2;

	% Tau evolution -2
	rho_L1=(U_left_tau*rho_L1); %.*T_relax;	
	rho_L2=(U_left_tau*rho_L2); %.*T_relax;
%}	

   	
	% t1 evolution: left half 
        for k=1:n-1
        rho_L1= U_left_t1*rho_L1; 					% tau evolution   
        rho_L2= U_left_t1*rho_L2; 					% tau evolution 
	end

 	% 180 H pulse +x for N decoupling
   	%rho_L=(L_pulse_180_H12x)*rho_L;		

	% t1 evolution: right half 
        for k=1:n-1
        rho_L1= U_left_t1*rho_L1; 					% tau evolution    
        rho_L2= U_left_t1*rho_L2; 					% tau evolution   
	end


 	% Quadrature detection: for cos:90 H pulse +x and 90 N pulse +x; for sin:H pulse +x and 90 N pulse +y
 	
   	rho_L1=(L_pulse_90_H1x * L_pulse_90_N1x)*rho_L1;
   	rho_L2=(L_pulse_90_H2x * L_pulse_90_N2x)*rho_L2;   	 	

	% Tau evolution -3
	rho_L1=(U_left_tau*rho_L1); %.*T_relax;
	rho_L2=(U_left_tau*rho_L2); %.*T_relax;
	
	
 	% 180 H pulse +x and 180 N pulse +x 
   	rho_L1=(L_pulse_180_H1x * L_pulse_180_N1x)*rho_L1;	 
   	rho_L2=(L_pulse_180_H2x * L_pulse_180_N2x)*rho_L2;
   	
   	
	% Tau evolution -4
	rho_L1=(U_left_tau*rho_L1); %.*T_relax;
	rho_L2=(U_left_tau*rho_L2); %.*T_relax;
	
	
	% FID acquisition	
        for k=1:td2							% t2
         fid_cos_1(n,k)=(conj(detop_L)*rho_L1);   		%fid(t1,t2)
         fid_cos_2(n,k)=(conj(detop_L)*rho_L2);   		%fid(t1,t2)

         fid_cos(n,k)=fid_cos_1(n,k)+fid_cos_2(n,k);   		%fid(t1,t2)
                 
	 rho_L1=(U_left_t2*rho_L1); %.*T_relax;
	 rho_L2=(U_left_t2*rho_L2); %.*T_relax;
	 
  	 fprintf(outfile_1,'%d\t%d\t%f\t%f\n',k,n,real(fid_cos(n,k)),imag(fid_cos(n,k)));
        end


% Sin modulation							% t1

    	rho_L1= sigma0_L;  
     	rho_L2= sigma0_L; 
     	
	% 90 H pulse +x
        rho_L1= L_pulse_90_H1x*rho_L1;			
        rho_L2= L_pulse_90_H2x*rho_L2;		

	% Tau evolution -1
	rho_L1=(U_left_tau*rho_L1); %.*T_relax;	
	rho_L2=(U_left_tau*rho_L2); %.*T_relax;	
   	
 	% 180 H pulse +x and 180 N pulse +x 
   	rho_L1=(L_pulse_180_H1x * L_pulse_180_N1x)*rho_L1;		
   	rho_L2=(L_pulse_180_H2x * L_pulse_180_N2x)*rho_L2;

	% Tau evolution -2
	rho_L1=(U_left_tau*rho_L1); %.*T_relax;	
	rho_L2=(U_left_tau*rho_L2); %.*T_relax;
		
 	% 90 H pulse +y and 90 N pulse +x 
   	rho_L1=(L_pulse_90_H1y * L_pulse_90_N1x)*rho_L1;	
   	rho_L2=(L_pulse_90_H2y * L_pulse_90_N2x)*rho_L2;

%{
		% CPMG RD element left half
		for cpmg_loop=1:N_CPMG_Pulse
			rho_L1=U_left_cpmg_tr_tau*rho_L1; 	% T_relax - tau;	
			rho_L2=U_left_cpmg_tr_tau*rho_L2; 	% T_relax - tau;			
			
		   	rho_L1=L_pulse_180_N1y*rho_L1; 	% pi pulse y
		   	rho_L2=L_pulse_180_N2y*rho_L2; 	% pi pulse y
		   	
		   	rho_L1=U_left_cpmg_tr_tau*rho_L1; 	% T_relax- tau;		   			   	
		   	rho_L2=U_left_cpmg_tr_tau*rho_L2; 	% T_relax- tau;	
		end

		
		%  CPMG : P element 
			rho_L1=U_left_cpmg_tau*rho_L1; 	% T_relax - tau;	
			rho_L2=U_left_cpmg_tau*rho_L2; 	% T_relax - tau;

		   	rho_L1=(L_pulse_90_H1y*L_pulse_180_H1x*L_pulse_90_H1y)*rho_L1; 	% p element pulses H		
		   	rho_L2=(L_pulse_90_H2y*L_pulse_180_H2x*L_pulse_90_H2y)*rho_L2; 	% p element pulses H	

		   	rho_L1=(L_pulse_180_N1x)*rho_L1; 	% p element pulses N			   		
		   	rho_L2=(L_pulse_180_N2x)*rho_L2; 	% p element pulses N	

			rho_L1=U_left_cpmg_tau*rho_L1; 	% T_relax - tau;		   		
			rho_L2=U_left_cpmg_tau*rho_L2; 	% T_relax - tau;	
			
		% CPMG RD element right half
		for cpmg_loop=1:N_CPMG_Pulse
			rho_L1=U_left_cpmg_tr_tau*rho_L1; 	% T_relax - tau;	
			rho_L2=U_left_cpmg_tr_tau*rho_L2; 	% T_relax - tau;	

		   	rho_L1=L_pulse_180_N1x*rho_L1; 	% pi pulse x
		   	rho_L2=L_pulse_180_N2x*rho_L2; 	% pi pulse x

		   	rho_L1=U_left_cpmg_tr_tau*rho_L1; 	% T_relax- tau;	
		   	rho_L2=U_left_cpmg_tr_tau*rho_L2; 	% T_relax- tau;	
		end


	% To recover IzSy from Sx

	% Tau evolution -1
	rho_L1=(U_left_tau*rho_L1); %.*T_relax;	
	rho_L2=(U_left_tau*rho_L2); %.*T_relax;	
   	
 	% 180 H pulse +x and 180 N pulse +x 
   	rho_L1=(L_pulse_180_H1x * L_pulse_180_N1x)*rho_L1;		
   	rho_L2=(L_pulse_180_H2x * L_pulse_180_N2x)*rho_L2;

	% Tau evolution -2
	rho_L1=(U_left_tau*rho_L1); %.*T_relax;	
	rho_L2=(U_left_tau*rho_L2); %.*T_relax;
%}	

   	
	% t1 evolution: left half 
        for k=1:n-1
        rho_L1= U_left_t1*rho_L1; 					% tau evolution   
        rho_L2= U_left_t1*rho_L2; 					% tau evolution 
	end

 	% 180 H pulse +x for N decoupling
   	%rho_L=(L_pulse_180_H12x)*rho_L;		

	% t1 evolution: right half 
        for k=1:n-1
        rho_L1= U_left_t1*rho_L1; 					% tau evolution    
        rho_L2= U_left_t1*rho_L2; 					% tau evolution   
	end


 	% Quadrature detection: for cos:90 H pulse +x and 90 N pulse +x; for sin:H pulse +x and 90 N pulse +y
 	
   	rho_L1=(L_pulse_90_H1x * L_pulse_90_N1y)*rho_L1;
     	rho_L2=(L_pulse_90_H2x * L_pulse_90_N2y)*rho_L2; 	 	

	% Tau evolution -3
	rho_L1=(U_left_tau*rho_L1); %.*T_relax;
	rho_L2=(U_left_tau*rho_L2); %.*T_relax;
		
 	% 180 H pulse +x and 180 N pulse +x 
   	rho_L1=(L_pulse_180_H1x * L_pulse_180_N1x)*rho_L1;	 
   	rho_L2=(L_pulse_180_H2x * L_pulse_180_N2x)*rho_L2;	   	
   	
	% Tau evolution -4
	rho_L1=(U_left_tau*rho_L1); %.*T_relax;
	rho_L2=(U_left_tau*rho_L2); %.*T_relax;
	
	% FID acquisition	
        for k=1:td2							% t2
         fid_sin_1(n,k)=(conj(detop_L)*rho_L1);   		%fid(t1,t2)
         fid_sin_2(n,k)=(conj(detop_L)*rho_L2);   		%fid(t1,t2)
         fid_sin(n,k)=fid_sin_1(n,k)+fid_sin_2(n,k);   		%fid(t1,t2)
                  
	 rho_L1=(U_left_t2*rho_L1); %.*T_relax;
	 rho_L2=(U_left_t2*rho_L2); %.*T_relax;
	 
  	 fprintf(outfile_1,'%d\t%d\t%f\t%f\n',k,n,real(fid_sin(n,k)),imag(fid_sin(n,k)));
        end


%pause

end % t1 loop

fclose(outfile_1);


File_str=sprintf('CPMG_RD_%dPulse_%.0f_kex_%.0fA-%.0fB_Pop',N_CPMG_Pulse,K,Pop_A,Pop_B);
fprintf(outfile_2,'\n# Expt: %d [TR(ms): %.1f pulses:%d pulses  Field:%0.3f Hz]',cpmg_expt,TR*1000,N_CPMG_Pulse,(1/(4*Tau_TR)));
fprintf(outfile_2,'\n# ----------------------------------------------------');
fprintf(outfile_2,'\ntxt2pipe.tcl < %s.asc > %s.fid -complex -xy -time -xN %d -yN %d -xT %d -yT %d -xMODE Complex -yMODE Complex -xSW %.3f -ySW %.3f -xOBS %.3f -yOBS %.3f -xCAR %.3f -yCAR %.3f -xLAB H -yLAB N -aq2D States -ndim 2',File_str,File_str,td2*2,td1*2,td2,td1,swh2,swh1,MHZ_1(1,1),abs(MHZ_1(1,2)),offset_ppm(1,1),offset_ppm(1,2));



fprintf(outfile_2,'\n\nnmrPipe  -in %s.fid \t\\',File_str);
fprintf(outfile_2,'\n| nmrPipe  -fn SP -off .9 -end .99 -pow 2 -c 1.0\t\\');
fprintf(outfile_2,'\n| nmrPipe  -fn ZF -size 2048\t\\');
fprintf(outfile_2,'\n| nmrPipe  -fn FT -auto\t\\');
fprintf(outfile_2,'\n| nmrPipe  -fn PS -p0 180 -p1 0.0 -di\t\\');
fprintf(outfile_2,'\n| nmrPipe  -fn TP\t\\');
fprintf(outfile_2,'\n| nmrPipe  -fn SP -off .9 -end .99 -pow 2 -c 1.0\t\\');
fprintf(outfile_2,'\n| nmrPipe  -fn ZF -size 2048\t\\');
fprintf(outfile_2,'\n| nmrPipe  -fn FT -auto\t\\');
fprintf(outfile_2,'\n| nmrPipe  -fn PS -p0 0 -p1 0 -di\t\\');
fprintf(outfile_2,'\n| nmrPipe  -fn REV\t\\');
fprintf(outfile_2,'\n| nmrPipe  -fn TP\t\\');
fprintf(outfile_2,'\n| nmrPipe  -out %s.ft2  -verb -ov',File_str);

fprintf(outfile_2,'\n\npipe2ucsf -12 %s.ft2 %s.ucsf\n\n',File_str,File_str);

% SHR phase processing

	% FOURIER TRANSFORM AND PLOTTING 
	% ------------------------------ 

	% States-TPPI based data processing
	%----------------------------------

	ft_1=(fft(fid_cos,[],2));
	ft_2=(fft(fid_sin,[],2));

	ft_3=real(ft_1)+complex(0,1)*real(ft_2);
	ft_4=fftshift(fft(ft_3,[],1));



	% Contour plot
	%--------------
	figure();
	zmin = min(min(real(ft_4)));
	zmin = 2.5*zmin;
	zmax =max(max(real(ft_4)));
	zinc = (zmax - zmin) / 100;
	zlevs = zmin:zinc:zmax;
	map=[1 0 0; 0 0 1];
	contour(f2,f1,real(ft_4),zlevs);

	set(gca,'CLim',[zmin zmax]);
	colormap(map);
	colorbar; 
	toc




end % cpmg_loop

fclose(outfile_2); % Pipe script







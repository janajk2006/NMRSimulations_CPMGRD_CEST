
%Calculates the steady state rho (density matrix) in Liouviille space, given Liouvillian, Sigmaeq=Sigma0

% transcribed by Jana.K
% Ackowledgement:  ScottA. Smith, Tilo Levante, T. S. Mahesh, Ilya Kuprov 
% =======================================================================

function [Sigma_SS]=Rho_steady_state(sigma0,Heff,L_RF);

% Steady state density matrix from RF relax Liovillian during pre-saturation period

%	L|s  > = R|s  >    -------> [L - L|S><E|]|s  > = R|s  > - |L1>
%         ss       eq		                   ss       eq

Dim_all=size(sigma0);
Dim=Dim_all(1,1); % symmetric matrices
DIM=Dim*Dim;



L_i=eye(Dim);
L_ROT=complex(0,1)*2*pi*(diprod(Heff,L_i)-diprod(L_i,Heff))+L_RF;


L_SS=zeros(DIM,1);


% RHS

sigma_eq=sigma0;

Emx_norm=eye(Dim,Dim)./(complex(Dim,0));
sigma_eq=sigma_eq+Emx_norm;
Sigmaeq_L=zeros(DIM,1);

%sigma_eq

% RHS
%			|s  > = R|s  > - |L1>
%			  eq       eq

R_SS=zeros(DIM,1); 	% Liouville superket of R_SS
L1_SS=zeros(DIM,1); 	% Liouville superket of R_SS


Sigmaeq_L=(reshape(transpose(sigma_eq),[],1));	
R_SS=L_RF*Sigmaeq_L;
L1_SS=L_ROT(:,1);
R_SS=R_SS-L1_SS;
R_SS_sub=R_SS(2:end,1);



Smx=zeros(Dim,Dim);
Smx(1,1)=1;
% Smx=Heff_V'*Smx(1,1)*Heff_V;

Emx=eye(Dim,Dim);
SMX=zeros(DIM,1); 
EMX=zeros(DIM,1); 

SMX=(reshape(transpose(Smx),[],1));
EMX=(reshape(transpose(Emx),[],1));


SE=zeros(DIM,DIM);
SE=SMX*transpose(EMX); %transpose(EMX)
X=L_ROT-L_ROT*SE;

% form a reduced X matrix as first column is zero
X=X(2:end,2:end); % (DIM-1),(DIM-1)

X_inv=inv(X);

L_SS_sub=X_inv*R_SS_sub;

%L_SS_sub
L_ss=zeros(Dim,Dim);

	m=1;
	for k=1:Dim
	for l=1:Dim
		if m~=1		% skip first element
		L_ss(k,l)=L_SS_sub(m-1,1);
		end
		m=m+1;
	end
	end

L_ss(1,1)=1-trace(L_ss);
L_ss=L_ss-Emx_norm;



Sigma_SS=(reshape(transpose(L_ss),[],1));



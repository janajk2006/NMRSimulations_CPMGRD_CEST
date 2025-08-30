clear;
clc;
pkg load statistics

% Number of spins
Nspins=1;

% Labels for the spins
Lspins=['I','S','X'];
Qspins=[1/2,1/2,1/2];

Lang_momentum=['i','x','y','z']; %  standard angular momentum operator labels

% loop over spins to generate the operators
outfile=fopen("pauli.m","w");

for i_1=1:Nspins



%Ie operator

	q=Qspins(1,i_1);
	Dim=floor(2*q+1);
	M=zeros(Dim,Dim);

	fprintf(outfile,"\n%ci=[",Lspins(1,i_1));
	for k=1:Dim
	M(k,k)=1;
	end % k
	
	for j=1:Dim	
	for k=1:Dim
	fprintf(outfile,"complex(%f,%f) ",real(M(j,k)),imag(M(j,k)));
	end % j
		if j<=Dim-1
		fprintf(outfile,"\n");
		else
		fprintf(outfile,"];\n");
		end
	end % k


% Iz operator

	q=Qspins(1,i_1);
	m=q;
	Dim=floor(2*q+1);
	M=zeros(Dim,Dim);

	fprintf(outfile,"\n%cz=[",Lspins(1,i_1));
	for k=1:Dim
	M(k,k)=m;
	m=m-1;
	end % k
	
	for j=1:Dim	
	for k=1:Dim
	fprintf(outfile,"complex(%f,%f) ",real(M(j,k)),imag(M(j,k)));
	end % j
		if j<=Dim-1
		fprintf(outfile,"\n");
		else
		fprintf(outfile,"];\n");
		end
	end % k


% Ix operator

	q=Qspins(1,i_1);
	m=q-1;
	q1=q*(q+1);
	Dim=floor(2*q+1);
	M=zeros(Dim,Dim);

	fprintf(outfile,"\n%cx=[",Lspins(1,i_1));
	for k=1:Dim-1
	M(k+1,k)=sqrt(q1-m*(m+1))/2;
	M(k,k+1)=sqrt(q1-m*(m+1))/2;
	m=m-1;
	end % k
	
	for j=1:Dim	
	for k=1:Dim
	fprintf(outfile,"complex(%f,%f) ",real(M(j,k)),imag(M(j,k)));
	end % j
		if j<=Dim-1
		fprintf(outfile,"\n");
		else
		fprintf(outfile,"];\n");
		end
	end % k

% Iy operator


	q=Qspins(1,i_1);
	m=q-1;
	q1=q*(q+1);
	Dim=floor(2*q+1);
	M=zeros(Dim,Dim);

	fprintf(outfile,"\n%cy=[",Lspins(1,i_1));
	for k=1:Dim-1
	M(k+1,k)=complex(0,sqrt(q1-m*(m+1))/2);
	M(k,k+1)=conj(complex(0,sqrt(q1-m*(m+1))/2));
	m=m-1;
	end % k
	
	for j=1:Dim	
	for k=1:Dim
	fprintf(outfile,"complex(%f,%f) ",real(M(j,k)),imag(M(j,k)));
	end % j
		if j<=Dim-1
		fprintf(outfile,"\n");
		else
		fprintf(outfile,"];\n");
		end
	end % k

% Ip operator

	q=Qspins(1,i_1);
	m=q-1;
	q1=q*(q+1);
	Dim=floor(2*q+1);
	M=zeros(Dim,Dim);

	fprintf(outfile,"\n%cp=[",Lspins(1,i_1));
	for k=1:Dim-1
	M(k,k+1)=sqrt(q1-m*(m+1));
	m=m-1;
	end % k
	
	for j=1:Dim	
	for k=1:Dim
	fprintf(outfile,"complex(%f,%f) ",real(M(j,k)),imag(M(j,k)));
	end % j
		if j<=Dim-1
		fprintf(outfile,"\n");
		else
		fprintf(outfile,"];\n");
		end
	end % k

% Im operator


	q=Qspins(1,i_1);
	m=q-1;
	q1=q*(q+1);
	Dim=floor(2*q+1);
	M=zeros(Dim,Dim);

	fprintf(outfile,"\n%cm=[",Lspins(1,i_1));
	for k=1:Dim-1
	M(k+1,k)=sqrt(q1-m*(m+1));
	m=m-1;
	end % k
	
	for j=1:Dim	
	for k=1:Dim
	fprintf(outfile,"complex(%f,%f) ",real(M(j,k)),imag(M(j,k)));
	end % j
		if j<=Dim-1
		fprintf(outfile,"\n");
		else
		fprintf(outfile,"];\n");
		end
	end % k


end %i_1 loop



fclose(outfile);



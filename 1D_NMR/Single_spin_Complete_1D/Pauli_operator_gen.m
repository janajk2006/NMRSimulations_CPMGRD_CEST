clear;
clc;
pkg load statistics

% Number of spins
Nspins=1;
Counter=power(power(2,Nspins),2);


% Labels for the spins
Lspins=['I','S','X'];

Lang_momentum=['i','x','y','z']; %  standard angular momentum operator labels

% loop over spins to generate the operators

Lang=[1:4];
for i=1:Nspins-1
Lang=[Lang 1:4];
end
A=nchoosek(Lang,Nspins); % change or repeat 1:4 equal to the number of spins
%A=nchoosek([1:4,1:4,1:4],Nspins); % change or repeat 1:4 equal to the number of spins
B=unique(A,'rows');
C=size(B);

outfile=fopen("nspins.m","w");

% cartesian operators
fprintf(outfile,"\npauli\n\n%% Cartesian operators: i,x,y,z\n\n");
for i=1:C(1,1)

	for j=1:Nspins
	fprintf(outfile,"%c%c",Lspins(1,j),Lang_momentum(1,B(i,j)));
	end
	fprintf(outfile,"=diprod(");

	for j=1:Nspins
	if j<=Nspins-1
	fprintf(outfile,"%c%c,",Lspins(1,j),Lang_momentum(1,B(i,j)));
	else
	fprintf(outfile,"%c%c);\n",Lspins(1,j),Lang_momentum(1,B(i,j)));
	end
	end

end

% raising, lowering operators
fprintf(outfile,"\n%% Ladder operators: i,+,-\n\n");


Lang_momentum=['i','p','m']; %  standard angular momentum ladder operator labels
Lang=[1:3];
for i=1:Nspins-1
Lang=[Lang 1:3];
end
A=nchoosek(Lang,Nspins); % change or repeat 1:4 equal to the number of spins
B=unique(A,'rows');
C=size(B);

for i=1:C(1,1)

	for j=1:Nspins
	fprintf(outfile,"%c%c",Lspins(1,j),Lang_momentum(1,B(i,j)));
	end
	fprintf(outfile,"=diprod(");

	for j=1:Nspins
	if j<=Nspins-1
	fprintf(outfile,"%c%c,",Lspins(1,j),Lang_momentum(1,B(i,j)));
	else
	fprintf(outfile,"%c%c);\n",Lspins(1,j),Lang_momentum(1,B(i,j)));
	end
	end

end




fclose(outfile);


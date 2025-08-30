function DP = diprod(varargin) 

% DP = diprod(M1,M2,...) 
% Gives direct product matrices M1,M2,... 
% T. S. Mahesh (20 May 2002) 


nomat=nargin;
DP=varargin{1}; 
for k=1:nomat-1 
	M=DP;
        N=varargin{k+1};
        [r s]=size(M); 
        [t u]=size(N); 
        for p=1:r 
        for q=1:s 
          DP((p-1)*t+1:(p-1)*t+t,(q-1)*u+1:(q-1)*u+u)=M(p,q)*N; 
        end 
        end 
end

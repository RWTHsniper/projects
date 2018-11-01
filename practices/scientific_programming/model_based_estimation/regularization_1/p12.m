clear;
A=[25.6 30.8; 23.2 32.6; 10.4 17.2];
[U,S,V]=svd(A);
y=[0;3;3];
St= S';
for i=1:2
St(i,i)=1/St(i,i);
end
x=V*St*U'*y;

%% TSVD
St2= St;
St2(2,2)=0;
T_x=V*St2*U'*y;

%% Tikhonov regularization

xa = zeros(size(x));
alpha = 5;
%alpha = 1*10^-5;
for i=1:size(xa)
   xa = xa + (S(i,i)^2/(S(i,i)^2+alpha^2)) * S(i,i)^-1 * V(:,i)* (U(:,i)'*y);
end


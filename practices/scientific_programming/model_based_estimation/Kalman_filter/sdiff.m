
%% This page shows how df/dx is done using Symbolic toolbox.

clear;
% symbolic differentiation
syms h1 h2 f1 f2 Fin T A_1 A_2 R_1 R_2;

f1 = h1 + T/A_1*(Fin - R_1*sqrt(h1));

 % df1/dh1
1 - (R_1*T)/(2*A_1*h1^(1/2))

% df1/dh2 
0

f2 = h2 + T/A_2*(R_1*sqrt(h1) - R_2*sqrt(h2));
diff(f2, h1)
(R_1*T)/(2*A_2*h1^(1/2))

diff(f2, h2)
1 - (R_2*T)/(2*A_2*h2^(1/2))




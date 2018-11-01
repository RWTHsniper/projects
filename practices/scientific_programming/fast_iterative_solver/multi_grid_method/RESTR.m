function [ u2h ] = RESTR(uh,Nc)
%RESTR Summary of this function goes here
%   Detailed explanation goes here

%I_h2H = zeros(Nc-1, 3+2*(Nc-1-1));
I_h2H = zeros(Nc-1, 2*Nc);

for i=1:Nc-1
   I_h2H(i,1+2*(i-1):3+2*(i-1)) = [1 2 1]; 
    
    
end

I_h2H = I_h2H/4;

u2h = I_h2H * uh;

end



clear;

N = 2^4;
Nc = N/2;
h = 1/N;

uh = zeros(N,N);

for i=1:N
    for j=1:N
x = i*h;
y = j*h;
        uh(i, j) = sin(2*pi*x)*sin(2*pi*y);
    end
end


%u2h = RESTR(uh,Nc);

I_h2H = zeros(Nc-1,2*Nc-1);

for i=1:Nc-1
   I_h2H(i,1+2*(i-1):3+2*(i-1)) = [1 2 1]; 
    
    
end

I_h2H = I_h2H/4;

u2h = I_h2H * uh(1:N-1,1:N-1);






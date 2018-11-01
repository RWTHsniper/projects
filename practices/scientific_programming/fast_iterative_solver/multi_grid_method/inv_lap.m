

%% inverse Laplacian operation
clear;

N = 2^2;
h= 1/N;
tol= 1e-10;
iteration = 0;

f= zeros(N+1,N+1);

u = zeros((N+1),(N+1)); % made it a little bit larger to accomodate outer boundary values
% boundary values are always zero
u_real = u;

for i=2:N+1
    for j=2:N+1
        x = (i-1)*h;
        y= (j-1)*h;
f(i,j) = 8*pi^2*sin(2*pi*x)*sin(2*pi*y);
u_real(i,j) = sin(2*pi*x)*sin(2*pi*y);
    end
end

u_v = zeros((N+1)*(N+1),1);
count = 1;
for i=1:(N+1)
    for j=1:(N+1)
u_v(count) = u_real(i,j);
count = count+1;
    end
end


Lap = zeros((N+1)*(N+1),(N+1)*(N+1));
Lap(1,1) = 4; 
Lap(1,2) = -1;
Lap (1, 1+(N+1)) = -1;
for i=2:(N+1)*(N+1) -1
   
    Lap(i,i) = 4;
       Lap(i, i+1) = -1;
       Lap(i, i-1) = -1;
    if (i+(N+1)) <=(N+1)*(N+1)
       Lap(i, i+(N+1)) = -1;
    end
    if (i-(N+1)) >= 1
       Lap(i, i-(N+1)) = -1;
    end

end

index = (N+1)*(N+1);

Lap(index,index) = 4;
Lap(index, index-1) = -1;
Lap(index, index - (N+1)) = -1;
Lap = Lap/h^2;

d_u_v = Lap* u_v;

u_v_2 = Lap\d_u_v;
% inverse laplacian works

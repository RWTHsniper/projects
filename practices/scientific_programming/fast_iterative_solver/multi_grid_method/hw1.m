
clear;

N = 2^5;
h= 1/N;
tol= 1e-10;
iteration = 0;

f= zeros(N+1,N+1);

u = zeros((N),(N)); % made it a little bit larger to accomodate outer boundary values
% boundary values are always zero
u_temp = u;
u_prev = u;
u_real = u;

for i=1:N-1
    for j=1:N-1
        x = i*h;
        y= j*h;
f(i,j) = 8*pi^2*sin(2*pi*x)*sin(2*pi*y);
u_real(i,j) = sin(2*pi*x)*sin(2*pi*y);
    end
end

%u = GS(u, f);

% Check result

u_temp = u_real - u;
norm = inf_norm(u_temp);
norm2 = inf_norm(u_real);
if (norm/norm2) < 1e-2
    disp('Solution converged');
else
    disp('Fucked up');
end
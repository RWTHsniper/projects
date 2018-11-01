
clear;

N = 2^7;
h= 1/N;
tol= 1e-10;
iteration = 0;

f= zeros(N+1,N+1);

u = zeros((N+1),(N+1)); % made it a little bit larger to accomodate outer boundary values
% boundary values are always zero
u_temp = u;
u_prev = u;
u_real = u;

for i=2:N
    for j=2:N
        x = (i-1)*h;
        y= (j-1)*h;
        u_real(i,j) = sin(2*pi*x)*sin(2*pi*y);
        f(i,j) = 8*pi^2*u_real(i,j);
    end
end

% Check result



% Laplacian operator done. u_temp = L * u_real.
% need to compare u_temp and f.
u = u_real;
for i=2:N
    for j=2:N
        
   u_temp(i,j)=-1/h^2*(u(i-1,j) - 4*u(i,j) + u(i+1,j) +u(i,j-1) + u(i,j+1));      
        
    end
end

u_temp2 = u_temp - f;

L = zeros((N+1)*(N+1), (N+1)*(N+1));

size = (N+1)*(N+1);
for i=1:size
L(i,i) = 4; 
end

for i=1:size-1
   L(i+1,i) = -1; 
end

for i=2:size
   L(i-1,i) = -1; 
end

for i=N+2:size
   L(i-(N+1),i) = -1;
   L(i,i-(N+1)) = -1;
end

L=L/h^2;

Lu = zeros(size,1);
c=1;
for i=1:(N+1)
   for j=1:(N+1)
      
       Lu(c,1)=u_temp(i,j);
       c=c+1;
   end
end

Lu_r = zeros(size,1);
c=1;
for i=1:(N+1)
   for j=1:(N+1)
      
       Lu_r(c,1)=u_real(i,j);
       c=c+1;
   end
end

%% L * Lu_r = f I guess.

disp(norm(u_temp-f)/norm(f));

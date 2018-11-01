function [ u ] = GS( u0, f )
%GS Summary of this function goes here
%   Detailed explanation goes here

iteration = 0;
u = u0;
tol = 1e-10;
N = length(u0);
h = 1/N;

while (1) 
iteration = iteration +1;
    u_prev = u;
    for j=1:N-1
   for i=1:N-1
       if (i > 1 && j >1)
      u(i,j) = 1/4*(f(i,j)*h^2 + u(i-1,j) + u(i,j-1) + u(i+1,j) + u(i, j+1)); 
       
       elseif ((i == 1) && j >1)
      u(i,j) = 1/4*(f(i,j)*h^2 +u(i,j-1) + u(i+1,j) + u(i, j+1)); 
       elseif (i > 1 && (j == 1))
      u(i,j) = 1/4*(f(i,j)*h^2 +u(i-1,j) + u(i+1,j) + u(i, j+1)); 
       end    
   end
    end

u_temp = u - u_prev;
norm = inf_norm(u_temp);
if norm < tol
    disp('Converged');
    disp(iteration);
    break;
end
end




end


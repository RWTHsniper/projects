function [ norm ] = inf_norm( A )
%INF_NORM Summary of this function goes here
%   Calculates infinite norm of a matrix

norm = 0;
for i=1:length(A)
    temp = 0;
   for j=1:length(A)
       temp = abs(A(i,j)) + temp;
   end
   
   if (norm < temp)
      norm = temp;
   end    
end

end


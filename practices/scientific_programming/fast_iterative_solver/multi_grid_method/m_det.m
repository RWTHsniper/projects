%% determinant
M=5;
A=rand(M,M);
tol = 1e-6;
A(1,1) = 0;
while abs(A(1,1)) < tol
    for i=2:M
    temp = A(1,:);
    A(1,:) = A(i,:);
    A(i,:) = temp;
  disp('exchanged');

      if A(1,1) > tol
        break;
      end
    end
end

for k=1:M-1
for i=k+1:M
   factor = A(i,k)/A(k,k);
    for j=k:M
        A(i,j) = A(i,j) - factor*A(k,j);
    end
    
    
end
end

mydet=1;
for i=1:M
   mydet = mydet * A(i,i); 
    
end

err = det(A) - mydet;


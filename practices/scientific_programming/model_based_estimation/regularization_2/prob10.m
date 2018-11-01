clear;
T=5;
L=pi;
n=100;
h=L/(n+1);
a=10^-3;
x=linspace(0,100,n)/100*L;

C=zeros(n,n);
for i=1:n
C(i,i)=2;
if i==n
   C(n-1,n)=-1;
   C(n,n-1)=-1;   
else
    C(i,i+1)=-1;
    C(i+1,i)=-1;
end
end

C=C*a/h^2;

[vect,val]=eig(C);
A=expm(-T*C);
u0=sin(x)+0.2*sin(8*x);
ut=expm(-T*C)*u0';
utmeas=ut+10^-2*randn(n,1);

plot(x,utmeas)
hold on
plot(x,ut,'*')
plot(x,u0)
legend('utmeas','ut','u0');
hold off
title('Initial value and Final one')

u0est=expm(T*C)*utmeas;

figure;
plot(x,u0est,'*')
hold on
plot(x,u0,'*')
legend('u0est','u0');
hold off
title('Estimation')
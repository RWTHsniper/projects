

t=zeros(6,1);
dt=0.15;
for i=1:6
    t(i) = dt* i;
end
X=[t exp(t) t.^3 sin(t)];
z= [1.2 0.6 1.6 0.9]';
y=X*z;
dy=10^(-2)*[-1 1 1 -0.5 -2 1]';
py=y+dy;

pz=X\py;

plot(z);
hold on;
plot(pz);
hold off;
title('zeta and p_zeta');

[U,S,V]=svd(X);
%[U2,S2,V2]=svd(A2);

tS=S';
r=rank(tS);
for i=1:r
   tS(i,i) = 1/S(i,i); 
end
tS(4,4)=0;
tz= V*tS*U'*py;

figure;
plot(z);
hold on;
plot(tz);
hold off;
title('zeta and t_zeta');

figure;
plot(y);
hold on;
plot(py);
plot(X*tz);
hold off;
title('data');
%% what happeens if sigma_reg(3) is zero? 
%% the estimation becomes less accurate.





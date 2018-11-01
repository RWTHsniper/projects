
clear;
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


%% normal equation
alpha=10^-3;
alpha=0.04218;
X2=(X'*X+alpha^2*eye(size(X'*X)));
y2= X'*y;
py2= X'*py;
pz2=X2\py2;
z2=X2\y2;

figure;
plot(pz2);
hold on;
plot(z);
hold off;
title('zeta and t zeta');

figure;
plot(y);
hold on;
plot(py);
plot(X*z2);
hold off;
title('datas');

alpha = linspace(10^-4,10^0,500);
z_err = zeros(size(alpha));
d_err = zeros(size(alpha));
for i=1:length(alpha)
X2=(X'*X+alpha(i)^2*eye(size(X'*X)));
y2= X'*y;
py2= X'*py;
pz2=X2\py2;
z2=X2\y2;

z_err(i)=norm(pz2 - z);
d_err(i)=norm(X*pz2 - py);
end

figure;
loglog(alpha, z_err);
hold on;
loglog(alpha, d_err);
hold off;
legend('z_err','d_err');
xlabel('alpha');




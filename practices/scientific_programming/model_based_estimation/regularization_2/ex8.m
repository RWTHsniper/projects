close all;
clear;
L = pi;
n= 100;
T = 5;
x=linspace(0,100,n)/100*L;
u0=sin(x)+0.2*sin(8*x);
alpha = 10^-3;
C = zeros(n,n);
for i=1:100
    C(i,i) = 2;
    if i < 100
        C(i,i+1) = -1;
        C(i+1, i) = -1;
    else
        C(i,i-1) = -1;
    end
end
h = L/(n-1);
C = C*alpha/h^2;
A = expm(-T*C);
ut=A*u0';
% applying homogeneous boundary condition
ut(1)=0;
ut(end)=0;
uT_meas = ut+10^-2*randn(n,1);
yp = uT_meas;

% building laplacian operator for a vector

%% (a)
beta = 1;
k = 250;
discrepancy = zeros(k,1);
flag = 0;
xi = zeros(n, 1);
for i=1:k
    xi = xi - beta* A'*(A*xi - yp);
    discrepancy(i) = norm(A*xi-yp);
    xi_v(i) = norm(xi);
    xi_cv(i) = norm(C*xi);
    if (discrepancy(i) <= norm(uT_meas - ut)) &&(flag ==0);
        xi_k = xi;
        kk = i;
        disp(kk);
        flag =1;
    end
end

dy = ones(k) * norm(uT_meas -ut);
if flag ==1
    fig1=figure( );
    semilogy(discrepancy, 'linewidth',5,'color','b');
    hold on
    semilogy(dy,'--','linewidth',5,'color','r');
    hold off
    grid;
    legend('Discrepancy curve','dy');
    xlabel('k');
    axis([0 100 0.05 0.5]);
    ylabel('Discrepancy');
    
    
    y=ut;
    yb = A* xi_k;
    fig2=figure();
    sizeofline=4;
    sizeofline2=6;
    plot(x, y,'linewidth',sizeofline);
    hold on;
    plot(x, yp,'*','linewidth',sizeofline2);
    plot(x, yb,'s','linewidth',sizeofline2);
    xlabel('x');
    ylabel('Final Temperature');
    grid;
    title('Final Temperature');
    legend('y (exact)','yp (measured)','yb (estimated)');
    
    fig3 = figure();
    plot(x, u0,'linewidth',sizeofline);
    hold on;
    plot(x, xi_k,'*','linewidth',sizeofline2);
    hold off
    legend('u0 (exact)','xi_k (estimated)');
    title('Initial Temperature');
end

%% (b)
starting_index=3;
figure;
loglog(discrepancy(starting_index:end),xi_v(starting_index:end),'*', 'linewidth',5);
title('Discrepancy vs xi-norm');
axis tight;
xlabel('Discrepancy');
ylabel('xi');

%% (c)


fig5=figure( );
loglog(discrepancy(starting_index:end,1),xi_cv(starting_index:end),'*', 'linewidth',5);
title('Discrepancy vs C-norm');
xlabel('Discrepancy');
ylabel('C-norm');
axis tight;


% Initial temperature and Final temperature estimation from (c)

% calculating xi when k* = 20
xi = zeros(n, 1);
k = 20;
for i=1:k
    xi = xi - beta* A'*(A*xi - yp);
    discrepancy(i) = norm(A*xi-yp);
    xi_v(i) = norm(xi);
    xi_cv(i) = norm(C*xi);
    if (discrepancy(i) <= norm(uT_meas - ut)) &&(flag ==0);
        xi_k = xi;
        kk = i;
        disp(kk);
        flag =1;
    end
end


y=ut;
yb = A* xi;
figure;
plot(x, y,'linewidth',sizeofline);
hold on;
plot(x, yp,'*','linewidth',sizeofline2);
plot(x, yb,'s','linewidth',sizeofline2);
xlabel('x');
ylabel('Final Temperature');
grid;
title('Final Temperature k^* from L-Curve');
legend('y (exact)','yp (measured)','yb (estimated)');


figure();
plot(x, u0,'linewidth',sizeofline);
hold on;
plot(x, xi_k,'*','linewidth',sizeofline2);
hold off
legend('u0 (exact)','xi_k (estimated)');
title('Initial Temperature');






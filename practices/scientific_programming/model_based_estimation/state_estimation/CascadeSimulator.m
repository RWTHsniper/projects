function [t,h1,h2]=CascadeSimulator(A,B,x0,sigma_x,sigma_y,tend,u,T)
%This function simulates continuously the tank cascade
%outputs: h1,h2 simulated measurements for the liquid height in the tanks

%State matrix
C1 = [1 0];
C2 = [0 1];

%System equations
ode = @(t,x) A * x + B*u(t) + sigma_x.* randn(size(x));

%Solve ode
%T sample time
[t,x] = ode45(ode,[0: T: tend],x0);


%Calculate the outputs
y_noise = sigma_y * randn(1,size(x,1));
h1 = C1 * x';
h2 = C2 * x' + y_noise;

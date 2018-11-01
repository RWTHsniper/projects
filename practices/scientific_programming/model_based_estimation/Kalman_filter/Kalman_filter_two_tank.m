% Delete all variables, close all figures,  clear the workspace
clear; close all; clc

%< ----------- System initial state and inputs -------------------->
A_1 = 2;                         % tank 1 cross sectional area, m^2
A_2 = 5;                         % tank 2 cross sectional area, m^2
tend = 40;                       % simulation end time
sigma_y = 0.1;                  % measurement noise std
sigma_x = 1e-2*[1,1]';                % process noise std
R_1 = 1;                         % resistance to the flow F_1
R_2 = 1;                         % resistance to the flow F_2
u_func   = @(t) sin(0.1*t)+1;    % Input function
T = 0.001;                        % sample time
t = 0:T:tend;                    % Time horizon for simulation
u = u_func(t);

x0 = [0 1]';                     % True Initial state of the process 

if 0
plot(t,u,'LineWidth',3);         % Plot input profile
set(gca,'FontSize',30)
grid on
xlabel('Time', 'FontSize', 24)
ylabel('Input', 'FontSize',24)
title('Input profile', 'FontSize',30)
end
%< ---------------------------------------------------------------->

%< -------------------- System matrices --------------------------->
% Discrete-time LTI matrices

F = zeros(2,1);
%< ---------------------------------------------------------------->

%< ----------------- Simulate the system -------------------------->
%[t,h1,h2]=CascadeSimulator(A,B,x0,sigma_x,sigma_y,tend,u_func,T);
[t,h1,h2]=CascadeSimulator(x0,sigma_x,sigma_y,tend,u_func,T);
yn = h2;
%< ---------------------------------------------------------------->

%< ------------------ Test Observability -------------------------->
%< ---------------------------------------------------------------->

%< ----------------- Kalman Filter -------------------------------->
x0_hat = [0.5 0.5]';                   % Estimate (guess) of initial state
x0_hat = [0.1 1.1]';                   % Estimate (guess) of initial state
%P0     = 1e-5*[1 0; 0 1];          % Covariance of initial state estimate
P0     = [1e-2 0; 0 1e-2];          % Covariance of initial state estimate
Q      = 1e-2*[1 0; 0 1];            % State noise covariance matrix
%Q      = 1e-1*[1 0; 0 1];            % State noise covariance matrix
R      = 1e2;                         % Measurement noise covariance matirx

xhat(:,1) = x0_hat;
P = P0;
Pout{1,1} = P0;

W=[1 0;0 1];
V=1;
Cd=[0 1];

for i=2:length(yn)
    % time update
    F(1) = sqrt(xhat(1,i-1))*R_1;
    F(2) = sqrt(xhat(2,i-1))*R_2;
    Ad = [(1 - (R_1*T)/(2*A_1*xhat(1,i-1)^(1/2))),  0
         (R_1*T)/(2*A_2*xhat(1,i-1)^(1/2)),  (1 -(R_2*T)/(2*A_2*xhat(2,i-1)^(1/2)))];
    
    xhat(2,i) = xhat(2,i-1) + T/A_2 * (F(1) - F(2));
    xhat(1,i) = xhat(1,i-1) + T/A_1 * (u(1,i-1) - F(1));
    
    P = Ad*P*Ad' + W*Q*W';
    
    % Measurement update
    K = P*Cd'*(Cd*P*Cd'+V*R*V')^-1;
    xhat(:,i) = xhat(:,i) + K * (yn(i) - xhat(2,i));
    P=(eye(size(P)) - K*Cd)*P*(eye(size(P)) - K*Cd)' + K*V*R*V'*K';
    
end
%< ---------------------------------------------------------------->

%< ------------------------ Plot the results ----------------------->
figure()
subplot(1,2,1)
plot(t,h1(1,:),'LineWidth', 3)
hold on
plot(t,xhat(1,:), '-.r','LineWidth',3)
grid on
set(gca,'FontSize',20)
xlabel('Time', 'FontSize',20)
ylabel('h1','FontSize',20)

subplot(1,2,2)
plot(t,h2(1,:),'LineWidth', 3)
hold on
plot(t,xhat(2,:), '-.r','LineWidth',3)
grid on
set(gca,'FontSize',20)
xlabel('Time','FontSize',20)
ylabel('h2','FontSize',20)

legend('Process', 'Kalman Estimate')
%< ----------------------------------------------------------------->
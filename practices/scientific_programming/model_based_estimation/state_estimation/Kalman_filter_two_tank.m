% Delete all variables, close all figures,  clear the workspace
clear; close all; clc

%< ----------- System initial state and inputs -------------------->
A_1 = 2;                         % tank 1 cross sectional area, m^2
A_2 = 5;                         % tank 2 cross sectional area, m^2
tend = 40;                       % simulation end time
sigma_y = 0.01;                  % measurement noise std
sigma_x = [1,1]';                % process noise std
R_1 = 1;                         % resistance to the flow F_1
R_2 = 1;                         % resistance to the flow F_2
u_func   = @(t) sin(0.1*t)+1;    % Input function
T = 0.01;                        % sample time
t = 0:T:tend;                    % Time horizon for simulation
u = u_func(t);

x0 = [0 1]';                     % True Initial state of the process 

plot(t,u,'LineWidth',3);         % Plot input profile
set(gca,'FontSize',30)
grid on
xlabel('Time', 'FontSize', 24)
ylabel('Input', 'FontSize',24)
title('Input profile', 'FontSize',30)

%< ---------------------------------------------------------------->

%< -------------------- System matrices --------------------------->
A = [-1/(A_1*R_1) 0;1/(A_2*R_1) -1/(A_2*R_2)];  % State matrix
B = [1/A_1; 0];                                 % Input matrix
C = [0 1];                                      % Output matrix
% Discrete-time LTI matrices
Ad = T*A + eye(size(A));                       % Discrete State matrix
Bd = T*B;                                      % Discrete Input matrix
Cd = C;
Dd = 0;
%< ---------------------------------------------------------------->

%< ----------------- Simulate the system -------------------------->
[t,h1,h2]=CascadeSimulator(A,B,x0,sigma_x,sigma_y,tend,u_func,T);
yn = h2;
%< ---------------------------------------------------------------->

%< ------------------ Test Observability -------------------------->
Od = [Cd Cd*Ad];
rank(Od)
%< ---------------------------------------------------------------->

%< ----------------- Kalman Filter -------------------------------->
x0_hat = [0.5 0.5]';                   % Estimate (guess) of initial state
P0     = [0.01 0; 0 0.01];          % Covariance of initial state estimate
Q      = [0.1 0; 0 0.1];            % State noise covariance matrix
R      = 1;                         % Measurement noise covariance matirx

xhat(:,1) = x0_hat;
P = P0;
Pout{1,1} = P0;

for i=2:length(yn)
    % time update
    xhat(:,i) = Ad*xhat(:,i-1) + Bd*u(:,i-1);
    P = Ad*P*Ad' + Q;
    
    % Measurement update
    K = P*Cd'*(Cd*P*Cd'+R)^-1;
    xhat(:,i) = xhat(:,i) + K * (yn(i) - Cd*xhat(:,i) -Dd*u(:,i));
    P=(eye(size(P)) - K*Cd)*P*(eye(size(P)) - K*Cd)' + K*R*K';
    
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
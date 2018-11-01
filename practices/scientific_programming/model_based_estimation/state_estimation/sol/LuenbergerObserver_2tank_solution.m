% Delete all variables, close all figures,  clear the workspace
clear; close all; clc


%< ----------- System initial state and inputs -------------------->
A_1 = 2;                         % tank 1 cross sectional area, m^2
A_2 = 5;                         % tank 2 cross sectional area, m^2
tend = 40;                       % simulation end time
sigma_y = 0.01;                  % measurement noise std
sigma_x = [1, 1]';               % process noise std
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
xlabel('Time (s)', 'FontSize', 24)
ylabel('F_{in}(m^3/s)', 'FontSize',24)
title('Input profile', 'FontSize',30)

%< ---------------------------------------------------------------->

%< -------------------- System matrices --------------------------->
A = [-1/(A_1*R_1) 0;1/(A_2*R_1) -1/(A_2*R_2)];  % State matrix
B = [1/A_1; 0];                                 % Input matrix
C = [0 1];                                      % Output matrix
% Discrete-time LTI matrices
Ad = (eye(size(A))+T*A);                       % Discrete State matrix
Bd = T*B;                                      % Discrete Input matrix
Cd = C;
Dd = 0;
%< ---------------------------------------------------------------->

%< -------------------- Fixing the observer gain ------------------>

L = place(Ad', Cd', [0.95 0.995])'

%< ----------------- Simulate the system -------------------------->
[t,h1,h2]=CascadeSimulator(A,B,x0,sigma_x,sigma_y,tend,u_func,T);
yn = h2;
%< ---------------------------------------------------------------->

%< --------------------- State estimation with noise -------------->
x0hatn = [0.5 0.5]';
xhatn(:,1) = x0hatn;
yhatn(:,1) = Cd*xhatn(:,1);
for i = 2:length(t)
xhatn(:,i) = Ad * xhatn(:,i-1)+Bd * u(1,i-1) + L *(yn(:,i-1) - yhatn(:,i-1));
yhatn(:,i) = Cd * xhatn(:,i);
end
%< ---------------------------------------------------------------->

%< ------------------------ Plot the results ----------------------->
figure()
subplot(1,2,1)
plot(t,h1(1,:),'LineWidth', 3)
hold on
plot(t,xhatn(1,:), '-.k','LineWidth',3)
grid on
set(gca,'FontSize',20)
xlabel('Time(s)', 'FontSize',20)
ylabel('h_1 (m)','FontSize',20)

subplot(1,2,2)
plot(t,h2(1,:),'LineWidth', 3)
hold on
plot(t,xhatn(2,:), '-.k','LineWidth',3)
grid on
set(gca,'FontSize',20)
xlabel('Time(s)','FontSize',20)
ylabel('h_2 (m)','FontSize',20)

legend('Process', 'State estimate')
%< ----------------------------------------------------------------->

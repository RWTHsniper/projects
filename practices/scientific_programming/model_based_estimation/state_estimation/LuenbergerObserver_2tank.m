% Delete all variables, close all figures,  clear the workspace
clear; close all; clc


%< ----------- System initial state and inputs -------------------->
A_1 = 2;                         % tank 1 cross sectional area, m^2
A_2 = 5;                         % tank 2 cross sectional area, m^2
tend = 40;                       % simulation end time
sigma_y = 0.01;                  % measurement noise std
sigma_x = [1,1]';            % process noise std
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
Ad = ;                       % Discrete State matrix
Bd = ;                                      % Discrete Input matrix
Cd = ;
Dd = ;
%< ---------------------------------------------------------------->

%< -------------------- Fixing the observer gain ------------------>

L = place( , , )'

%< ----------------- Simulate the system -------------------------->
[t,h1,h2]=CascadeSimulator(,,,,,,,);
yn = h2;
%< ---------------------------------------------------------------->

%< --------------------- State estimation ------------------------->
x0hatn = [0.5 0.5]';

%< ---------------------------------------------------------------->

%< ------------------------ Plot the results ----------------------->
figure()
subplot(1,2,1)
plot(t,,'LineWidth', 3)
hold on
plot(t,, '-.k','LineWidth',3)
grid on
set(gca,'FontSize',20)
xlabel('', 'FontSize',20)
ylabel('','FontSize',20)

subplot(1,2,2)
plot(t,,'LineWidth', 3)
hold on
plot(t,, '-.k','LineWidth',3)
grid on
set(gca,'FontSize',20)
xlabel('','FontSize',20)
ylabel('','FontSize',20)

legend('Process', 'State estimate')
%< ----------------------------------------------------------------->
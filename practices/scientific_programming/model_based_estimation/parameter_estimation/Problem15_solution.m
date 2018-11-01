function [p_opt feval] = Problem15_solution()
% nonlinear parameter estimation
% using fmincon

parameter_estimation_data;                                                  % load the experimental data

figure
plot(u, y_tilde, 'LineWidth', 2)
xlabel('u')
ylabel('y\_tilde')
xlim([0 3])

options = optimset( 'Algorithm', 'active-set');                             % specify the optimization algorithms

p0 = [0;0];                                                                 % initial guess for the optimization
[p_opt feval] = fmincon(@(p)myfun(p,u, y_tilde),p0,[],[],[],[],...            % Call fmincon
    [0 ;0],[10;10],[], options);
%% Part iii: With initial guess p0 = [10;10]; Different solution is obtained

end


% function to calculate the objective function
function sse = myfun(p, u, y_tilde)

y = size(y_tilde);                                                         % specifying the size of z2_pred
sse = 0;
for i = 1:length(u)
    y(i) = p(1) + 1/(u(i) - p(2));                                         % calculate the predicted measurement
    sse = sse + (y_tilde(i) - y(i))^2;                                     % calculate the sum of squared errors
end

end
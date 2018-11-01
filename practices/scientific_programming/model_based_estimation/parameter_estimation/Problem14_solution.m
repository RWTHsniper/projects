clear; close all; clc
% number of Parameters
numParam = 3;

% Inputs
% Use -->u = startValue : stepSize: endValue;<--
% u = 
u = [0:0.1:1]';

% The measurements:
y_n = [1.06458112421035,0.893487378302637,1.09852277251498,1.17183112728074,1.03949589535870,0.974481255458338,0.823057849447003,0.918651740367926,0.685508877340526,0.566571755467876,0.569930393861726]';

std = [0.0645811242103507,0.146512621697363,0.0385227725149822,0.111831127280739,0.000504104641302350,0.0255187445416618,0.116942150552997,0.0586517403679255,0.0744911226594734,0.0734282445321243,0.0699303938617260];

%% Calculation of the X matrix in the equation X.theta = Y

X = zeros(length(u), numParam );
for i = 1:length(u)
    X(i,:) = [1 u(i) u(i)^2];
end

W = 1./std;

XW = diag(W) * X;
YW = diag(W) * y_n;

thetaHatW = (XW' * XW)\XW' * y_n
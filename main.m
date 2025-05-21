%% Me 469 Project
% Spring 2025
close all; clear; clc;

%% Params
% Body Params
lt = 0.47;
lh = 0.21;
mt = 60.5;
mh = 4.3;

% Spring/Damper Params
kb = 6126;
kn = 0;
bb = 1200;
bn = 0;

% Other Variables
g = 9.81; % m/s
l_body = lt*lh;


%% Torso Simulation
% Equation Break-up
A_t_12 = ((3/2) * (mt*lt*g - 2*kb)) / (mt*lt^2);
A_t_22 = - (3*bb) / (mt*lt^2);

% State Space
A_t = [0 1;
     A_t_12 A_t_22];
B_t = [0; 0];
C = eye(2,2);
D = [0; 0];

IC = [50; 0];

% Simulation
sys = ss(A_t,B_t,C,D);
t = 0:0.001:2;
n = length(t);
u = zeros(1, n);
x_t = lsim(sys, u, t, IC);

theta_t = x_t(:,1);
theta_t_dot = x_t(:,2);
theta_t_ddot = A_t_12 * theta_t + A_t_22 * theta_t_dot; % CHECK THIS



% %% Animation of Torso with Spring and Damper
% animateTorso(t, x_t, lt)
% 
% % Final plot showing the trajectory
% figure;
% plot(t, theta_t, 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('Angle (degrees)');
% title('Torso Angle vs Time');
% grid on;


%% Head Simulation
x_h = zeros(n, 2);
x_h(1, :) = [0; 0];

% Conversion
theta_t_rad = theta_t * pi/180;
theta_t_dot_rad = theta_t_dot * pi/180;
theta_t_ddot_rad = theta_t_ddot * pi/180;

% Numerical Integration
for i = 1:n-1
     % Current state
     theta_h = x_h(i,1) * pi/180;
     theta_h_dot = x_h(i,2) * pi/180;
    
     theta_t_i = theta_t_rad(i);
     theta_t_dot_i = theta_t_dot_rad(i);
     theta_t_ddot_i = theta_t_ddot_rad(i);
    
     % Calculate coefficients for this time step
     den = mh*lh^2;
     A_h_12 = - (kn + mh*l_body*(theta_t_ddot_i*sin(theta_t_i)+theta_t_dot_i^2*cos(theta_t_i))-mh*g*lh) / den;
     A_h_22 = - bn / den;
     A_h = [0, 1; 
          A_h_12, A_h_22];
    
    % Calculate forcing term (input from torso motion)
     B_h = [0;
          lt / lh];
    
     u = theta_t_ddot_i*sin(theta_t_i)+theta_t_dot_i^2*cos(theta_t_i);
     % Current state derivative: ẋ = Ax + B
     x_h_dot = A_h * [theta_h; theta_h_dot] + B_h * u;
    
     % Numerical integration (Euler method)
     dt = t(i+1) - t(i);
     x_h(i+1,1) = x_h(i,1) + x_h_dot(1) * dt;
     x_h(i+1,2) = x_h(i,2) + x_h_dot(2) * dt;
end

% Extract results
theta_h = x_h(:,1) * (180/pi);         % Head angle (degrees)
theta_h_dot = x_h(:,2)  * (180/pi);     % Head angular velocity (degrees/s)

%% Animation of Torso and Head
% Animation of both Torso and Head
figure;
animateTorsoAndHead(t, theta_t, theta_h, lt, lh)

% Comparison plot showing both angles
figure;
plot(t, theta_t, 'LineWidth', 2); hold on;
plot(t, theta_h, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Angle (degrees)');
title('Torso and Head Angles');
legend('Torso', 'Head');
grid on;
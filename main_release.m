%% ME 469 Project – Rear-End Collision Whiplash Model with Animation
close all; clear; clc;

%% Parameters
% Body parameters
lt = 0.47;    % Torso length (m)
lh = 0.21;    % Head/neck length (m)
mt = 60.5;    % Torso mass (kg)
mh = 4.3;     % Head mass (kg)

% Neck stiffness & damping
kb = 6126;         % Torso stiffness
kn_initial = 0.5;  % Very low pre-impact neck stiffness
bn_initial = 1;  % Very low pre-impact neck damping
kn_final = 35;     % High post-impact stiffness
bn_final = 0.5;      % High post-impact damping

% Other
g = 9.81;
l_body = lt * lh;

%% Torso Simulation
A_t_12 = ((3/2) * (mt * lt * g - 2 * kb)) / (mt * lt^2);
A_t_22 = - (3 * 1200) / (mt * lt^2);

A_t = [0 1; A_t_12 A_t_22];
B_t = [0; 0];
C = eye(2);
D = [0; 0];
IC = [50; 0];  % Initial torso angle in degrees

% Time setup
t = 0:0.001:2;
n = length(t);
u = zeros(1, n);

% Run simulation
sys = ss(A_t, B_t, C, D);
x_t = lsim(sys, u, t, IC);

theta_t = x_t(:,1);           % degrees
theta_t_dot = x_t(:,2);       % deg/s
theta_t_ddot = A_t_12 * theta_t + A_t_22 * theta_t_dot;

% Convert to radians
theta_t_rad = theta_t * pi/180;
theta_t_dot_rad = theta_t_dot * pi/180;
theta_t_ddot_rad = theta_t_ddot * pi/180;

%% Head Simulation (Euler Integration)
x_h = zeros(n, 2);
x_h(1, :) = [5 * pi/180, 0];  % Start at 5 deg forward tilt

head_released = false;
kn_active = kn_initial;
bn_active = bn_initial;

for i = 1:n-1
    % Current state
    theta_h = x_h(i,1);
    theta_h_dot = x_h(i,2);

    % Torso state
    theta_t_i = theta_t_rad(i);
    theta_t_ddot_i = theta_t_ddot_rad(i);

    % Trigger release when torso is about 45 deg
    if ~head_released && theta_t(i) < 45
        head_released = true;
        kn_active = kn_final;
        bn_active = bn_final;
    end

    % System dynamics
    den = mh * lh^2;
    A_h_12 = - (kn_active + mh * l_body * (theta_t_ddot_i * sin(theta_t_i)) - mh * g * lh) / den;
    A_h_22 = - bn_active / den;
    A_h = [0, 1; A_h_12, A_h_22];

    % Forcing term (only angular acceleration)
    B_h = [0; lt / lh];
    if ~head_released
        u = -theta_t_ddot_i * sin(theta_t_i);  % reverse input pre-impact
    else
        u = theta_t_ddot_i * sin(theta_t_i);   % normal input post-impact
    end

    % Euler integration
    x_h_dot = A_h * [theta_h; theta_h_dot] + B_h * u;
    dt = t(i+1) - t(i);
    x_h(i+1,1) = x_h(i,1) + x_h_dot(1) * dt;
    x_h(i+1,2) = x_h(i,2) + x_h_dot(2) * dt;
end

% Convert to degrees
theta_h = x_h(:,1) * 180/pi;
theta_h_dot = x_h(:,2) * 180/pi;

%% Animation (uses your existing function)
figure;
animateTorsoAndHead(t, theta_t, theta_h, lt, lh)

%% Plot
figure;
plot(t, theta_t, 'b', 'LineWidth', 2); hold on;
plot(t, theta_h, 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Angle (degrees)');
title('Torso and Head Angles – Rear-End Collision');
legend('Torso', 'Head');
grid on;
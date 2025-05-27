%% ME 469 Project – Rear-End Collision Whiplash Model with Animation
close all; clear; clc;

%% Parameters
% Body parameters
lt = 0.47;    % Torso length (m)
lh = 0.21;    % Head/neck length (m)
mt = 60.5;    % Torso mass (kg)
mh = 4.3;     % Head mass (kg)
nfig = 0;

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

    % Store the derivatives for this time step
    x_h_dot_history(i,:) = x_h_dot';

    dt = t(i+1) - t(i);
    x_h(i+1,1) = x_h(i,1) + x_h_dot(1) * dt;
    x_h(i+1,2) = x_h(i,2) + x_h_dot(2) * dt;
end

% Calculate the final time step derivative (can use previous step or set to zero)
x_h_dot_history(n,:) = x_h_dot_history(n-1,:);  % Use last calculated value for final time step


% Convert to degrees
theta_h_deg = x_h(:,1) * 180/pi;
theta_h_dot_deg = x_h(:,2) * 180/pi;
theta_h_ddot_deg = x_h_dot_history(:,2) * 180/pi;  % Angular acceleration in deg/s²
theta_h_ddot = theta_h_ddot_deg * pi/180;

g_force = (lh * theta_h_ddot)  / 9.8;



%% Animation (uses your existing function)
nfig = nfig + 1;
animateTorsoAndHead(t, theta_t, theta_h_deg, lt, lh)

%% Plot angles, velocities and accelerations
nfig = nfig + 1;
figure(nfig);
subplot(4,1,1);
plot(t, theta_t, 'b', 'LineWidth', 1); 
hold on;
plot(t, theta_h_deg, 'r', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Angular Acceleration (deg/s²)', 'FontSize', 12);
title('Torso and Head Angles');
legend('Torso', 'Head', 'best');
grid on;

subplot(4,1,2);
plot(t, theta_t_dot, 'b', 'LineWidth', 1); 
hold on;
plot(t, theta_h_dot_deg, 'r', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Angular Acceleration (deg/s²)', 'FontSize', 12);
title('Torso and Head Angular Velocities');
grid on;

subplot(4,1,3);
plot(t, theta_t_ddot, 'b', 'LineWidth', 1); 
hold on;
plot(t, theta_h_ddot, 'r', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Angular Acceleration (deg/s²)', 'FontSize', 12);
title('Torso and Head Angular Accelerations');
grid on;

subplot(4,1,4);
plot(t, g_force, 'r', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 12);
ylabel('G-Force (N)', 'FontSize', 12);
title('Head G-Force Experienced');
legend('Head')
grid on;

% Adjust figure size
set(gcf, 'Position', [100 100 800 800]);

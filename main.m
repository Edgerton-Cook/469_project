1%% Me 469 Project
% Spring 2025
close all; clear; clc;
nfig = 0;
set(0, 'DefaultFigureWindowStyle', 'docked')

%% Params
% Body Params
lt = 0.47;
lh = 0.21;
mt = 60.5;
mh = 4.3;

% Spring/Damper Params
kb = 7000;
kn = 100;
bb = 150;
bn = 10;

% Other Variables
g = 9.81; % m/s
l_body = lt*lh;

% Booleans
torso_animation = 0;
full_animation = 0;



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
t = 0:0.001:0.5;
n = length(t);
u_vec = zeros(1, n);
x_t = lsim(sys, u_vec, t, IC);

theta_t = x_t(:,1);
theta_t_dot = x_t(:,2);
theta_t_ddot = A_t_12 * theta_t + A_t_22 * theta_t_dot; % CHECK THIS


if torso_animation
    % Animation of Torso with Spring and Damper
    animateTorso(t, x_t, lt)
    
    % Final plot showing the trajectory
    figure;
    plot(t, theta_t, 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Angle (degrees)');
    title('Torso Angle vs Time');
    grid on;
end


%% Head Simulation
x_h = zeros(n, 2);
x_h(1, :) = [0; 0] * (pi/180);

% Conversion
theta_t_rad = theta_t * pi/180;
theta_t_dot_rad = theta_t_dot * pi/180;
theta_t_ddot_rad = theta_t_ddot * pi/180;

% Numerical Integration
for i = 1:n-1
     % Current state
     theta_h = x_h(i,1);
     theta_h_dot = x_h(i,2);
    
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

%% Plot angles, velocities, accelerations, and jerk
nfig = nfig + 1;
figure(nfig);
subplot(2,1,1);
plot(t, theta_t, 'b', 'LineWidth', 1); 
hold on;
plot(t, theta_h_deg, 'r', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Angle (deg)', 'FontSize', 12);
title('Torso and Head Angles');
legend('Torso', 'Head', 'best');
grid on;

subplot(2,1,2);
plot(t, g_force, 'r', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 12);
ylabel('G-Force', 'FontSize', 12);
title('Head G-Force Experienced');
grid on;

% Adjust figure size
set(gcf, 'Position', [100, 100, 900, 1000]);

if full_animation
    % Animation of Head with Torso
    figure;
    animateTorsoAndHead(t, theta_t, theta_h_deg, lt, lh);
end

clearvars -except x_h theta_h_ddot t nfig mh lh u_vec g_force
%% Second Simulation

% Find initial condition when head changes direction
found_change = false;
for i = 2:length(t)
    if g_force(i) > 16  % Sign change
        v_i = x_h(i,2);  % Angular velocity at this point
        t_crit = t(i);
        break
    end
end

IC = [0; v_i];  % Initial condition [angle; velocity]


% Define expanded parameter values to test
kr_values = [50, 100, 200, 400, 800, 1600];  % Extended range of spring constants
br_values = [1, 5, 10, 20, 40, 80];          % Extended range of damping values

% Create figure with 2 subplots
figure('Position', [100, 100, 1200, 800]);

% Subplot for position (angle)
subplot(2,1,1);
hold on;
title('Head Angle During Recovery', 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Angle (degrees)', 'FontSize', 12);
grid on;

% Subplot for g-force
subplot(2,1,2);
hold on;
title('G-Force During Recovery', 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 12);
ylabel('G-Force', 'FontSize', 12);
grid on;

% Store max g-force for each parameter set
max_g_force = zeros(length(kr_values), length(br_values));
legend_text = cell(length(kr_values) * length(br_values), 1);
plot_idx = 1;

% Define base colors for each kr value (distinct hues)
kr_base_colors = [
    1, 0.4, 0.4;     % Light Red for kr=50
    1, 0, 0;         % Red for kr=100
    0, 0, 1;         % Blue for kr=200
    0, 0.7, 0;       % Green for kr=400
    0.8, 0, 0.8;     % Purple for kr=800
    0.5, 0.5, 0;     % Olive for kr=1600
];

% Run simulations for each parameter combination
for i = 1:length(kr_values)
    kr = kr_values(i);
    
    % Get base color for this kr value
    base_color = kr_base_colors(i,:);
    
    for j = 1:length(br_values)
        br = br_values(j);
        
        % Create state-space model
        A = [0, 1; -kr/mh, -br/mh];
        B = [0; 1/mh];
        C = eye(2);
        D = [0; 0];
        
        sys = ss(A, B, C, D);
        
        % Simulate recovery
        [y, ~, x] = lsim(sys, u_vec, t, IC);
        
        % Extract states
        theta = y(:,1);      % Angle in radians
        theta_dot = y(:,2);  % Angular velocity in rad/s
        
        % Convert angle to degrees for plotting
        theta_deg = theta * 180/pi;
        
        % Calculate acceleration and g-force
        theta_ddot = -kr/mh * theta - br/mh * theta_dot;
        g_force = abs(lh * theta_ddot) / 9.81;
        max_g_force(i,j) = max(g_force);
        
        % Calculate brightness based on br value
        % Lower br = darker, Higher br = lighter
        brightness = 0.3 + 0.7 * (j-1)/(length(br_values)-1);
        
        % Adjust the color based on the br value
        plot_color = base_color * brightness;
        
        % Plot angle in top subplot
        subplot(2,1,1);
        plot(t, theta_deg, 'LineWidth', 1, 'Color', plot_color);
        
        % Plot g-force in bottom subplot
        subplot(2,1,2);
        plot(t, g_force, 'LineWidth', 1, 'Color', plot_color);
        
        % Prepare legend entry
        legend_text{plot_idx} = sprintf('kr=%d, br=%d', kr, br);
        plot_idx = plot_idx + 1;
    end
end

% Create simplified legend entries grouped by kr
kr_legend_entries = cell(length(kr_values), 1);
for i = 1:length(kr_values)
    kr_legend_entries{i} = sprintf('kr=%d', kr_values(i));
end

% Add color legend to the first subplot
subplot(2,1,1);
h_legend = legend(kr_legend_entries, 'Location', 'eastoutside', 'FontSize', 10);
title(h_legend, 'Spring Constants');

% Add a small note about brightness
annotation('textbox', [0.85, 0.48, 0.15, 0.1], ...
    'String', {'Brightness indicates', 'damping (br):', 'Darker = lower br', 'Lighter = higher br'}, ...
    'EdgeColor', 'none', 'FontSize', 9);



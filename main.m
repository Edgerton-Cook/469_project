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
kn = 35500;
bb = 1200;
bn = 421;

% Other Variables
g = 9.81; % m/s
l_body = lt*lh;

% Plotting Colors
Colors = {
    [0.0000 0.4470 0.7410],
    [0.8500 0.3250 0.0980],
    [0.9290 0.6940 0.1250],
    [0.4940 0.1840 0.5560],
    [0.4660 0.6740 0.1880],
    [0.3010 0.7450 0.9330],
    [0.6350 0.0780 0.1840]
};

%% Torso Simulation
% Equation Break-up
den = mt * lt^2;
A_12 = (3/2) * (mt*lt*g - 2*kb) / den;
A_22 = - 3*bb / (mt*lt^2);

% State Space
A = [0 1;
     A_12 A_22];
B = [0; 0];
C = eye(2,2);
D = [0; 0];

IC = [50; 0];

% Simulation
sys = ss(A,B,C,D);
t = 0:0.001:2;
u = zeros(1, length(t));
y = lsim(sys, u, t, IC);

%% Animation of Torso with Spring and Damper
% Create figure
figure('Position', [100, 100, 800, 600]);
h_ax = axes('XLim', [-1.2*lt, 1.2*lt], 'YLim', [-0.2, 1.2*lt]);
axis equal;
grid on;
hold on;

% Initialize graphics objects
torso = line([0, 0], [0, 0], 'LineWidth', 1, 'Color', Colors{1});
spring = line([0, 0], [0, 0], 'LineWidth', 1, 'Color', Colors{2}, 'LineStyle', ':');
damper = line([0, 0], [0, 0], 'LineWidth', 1, 'Color', Colors{5}, 'LineStyle', '--');
base = line([-0.2, 0.2], [0, 0], 'LineWidth', 1, 'Color', 'k');
title('Torso Simulation with Spring and Damper');
xlabel('X position (m)');
ylabel('Y position (m)');

% Draw legend
legend('Torso', 'Spring', 'Damper', 'Base', 'Location', 'best');

% Animation loop
for i = 1:10:length(t)
    % Get current angle
    theta = y(i, 1) * pi/180;  % Convert to radians if your theta is in degrees
    
    % Calculate torso endpoint
    x_torso = lt * sin(theta);
    y_torso = lt * cos(theta);
    
    % Update torso position
    set(torso, 'XData', [0, x_torso], 'YData', [0, y_torso]);
    
    % Calculate spring attachment point (at top of torso)
    spring_x = [x_torso, 0];   % Spring is attached to a fixed point at x=0
    spring_y = [y_torso, lt];  % Spring is attached to a fixed point at y=lt
    
    % Update spring position
    set(spring, 'XData', spring_x, 'YData', spring_y);
    
    % Update damper position (parallel to spring)
    set(damper, 'XData', spring_x, 'YData', spring_y);
    
    % Add text for current time and angle
    delete(findobj(h_ax, 'Type', 'text'));
    text(-1.1*lt, 1.1*lt, sprintf('Time: %.2f s', t(i)), 'FontSize', 10);
    text(-1.1*lt, 1.0*lt, sprintf('Angle: %.2f deg', y(i,1)), 'FontSize', 10);
    
    drawnow;
    pause(0.01);
end

% Final plot showing the trajectory
figure;
plot(t, y(:,1), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Angle (degrees)');
title('Torso Angle vs Time');
grid on;

function animateTorso(t, y, lt)
%ANIMATETORSO Animates a torso with spring and damper
%   Displays an animation of a rod (torso) rotating around a fixed point
%   with spring and damper attachments
%
%   Input:
%       t  - time vector
%       y  - simulation output (column 1: angle in degrees, column 2: angular velocity)
%       lt - length of the torso

    % Create figure
    figure('Position', [100, 100, 800, 600]);
    h_ax = axes('XLim', [-1.2*lt, 1.2*lt], 'YLim', [-0.2, 1.2*lt]);
    axis equal;
    grid on;
    hold on;
    
    % Define colors for visualization
    Colors = {
        [0.0000 0.4470 0.7410], % Blue
        [0.8500 0.3250 0.0980], % Red
        [0.9290 0.6940 0.1250], % Yellow
        [0.4940 0.1840 0.5560], % Purple
        [0.4660 0.6740 0.1880], % Green
        [0.3010 0.7450 0.9330], % Light Blue
        [0.6350 0.0780 0.1840]  % Burgundy
    };

    % Initialize graphics objects
    torso = line([0, 0], [0, 0], 'LineWidth', 1, 'Color', Colors{1});
    spring = line([0, 0], [0, 0], 'LineWidth', 1, 'Color', Colors{5}, 'LineStyle', '-');
    damper = line([0, 0], [0, 0], 'LineWidth', 1, 'Color', Colors{2}, 'LineStyle', '-');
    base = line([-0.2, 0.2], [0, 0], 'LineWidth', 1, 'Color', 'k');
    
    % Add wall on left side with hatching
    wall_x = -1*lt;
    wall = line([wall_x, wall_x], [0, lt], 'LineWidth', 1, 'Color', 'k');
    % Add hatching
    for i = 0:0.05:lt
        line([wall_x, wall_x-0.05], [i, i+0.05], 'LineWidth', 0.5, 'Color', 'k');
    end
    
    % Add connection points on wall
    spring_wall_y = 0.99*lt;
    damper_wall_y = 0.98*lt;
    
    title('Torso Simulation with Spring and Damper', 'FontWeight', 'normal');
    xlabel('Position (m)');
    ylabel('Height (m)');
    
    % Draw legend
    legend('Torso', 'Spring', 'Damper', 'Base', 'Location', 'best');
    
    % Animation loop
    for i = 1:10:length(t)
        % Get current angle
        theta = y(i, 1) * pi/180;  % Convert to radians
        
        % Calculate torso endpoint
        x_torso = lt * sin(theta);
        y_torso = lt * cos(theta);
        
        % Update torso position
        set(torso, 'XData', [0, x_torso], 'YData', [0, y_torso]);
        
        % Calculate spring attachment point (on torso, 2/3 of the way up)
        spring_torso_x = 0.99 * x_torso;
        spring_torso_y = 0.99 * y_torso;
        
        % Calculate damper attachment point (on torso, 1/3 of the way up)
        damper_torso_x = 0.98 * x_torso;
        damper_torso_y = 0.98 * y_torso;
        
        % Update spring position
        set(spring, 'XData', [wall_x, spring_torso_x], 'YData', [spring_wall_y, spring_torso_y]);
        
        % Update damper position
        set(damper, 'XData', [wall_x, damper_torso_x], 'YData', [damper_wall_y, damper_torso_y]);
        
        % Add mass at end of torso
        % delete(findobj(h_ax, 'Tag', 'mass'));
        % mass_radius = 0.06;
        % th_circle = linspace(0, 2*pi, 30);
        % x_circle = x_torso + mass_radius * cos(th_circle);
        % y_circle = y_torso + mass_radius * sin(th_circle);
        % plot(x_circle, y_circle, 'k-', 'LineWidth', 1, 'Tag', 'mass');
        
        % Draw angle indicator
        % delete(findobj(h_ax, 'Tag', 'angle'));
        % r = 0.15;
        % th = linspace(0, theta, 20);
        % x_arc = r * sin(th);
        % y_arc = r * cos(th);
        % plot(x_arc, y_arc, 'k-', 'LineWidth', 1, 'Tag', 'angle');
        
        % Add text for current time and angle
        delete(findobj(h_ax, 'Tag', 'info_text'));
        text(-1.1*lt, 1.1*lt, sprintf('Time: %.2f s', t(i)), 'FontSize', 10, 'Tag', 'info_text');
        text(-1.1*lt, 1.0*lt, sprintf('Angle: %.2f deg', y(i,1)), 'FontSize', 10, 'Tag', 'info_text');
        
        drawnow;
        pause(0.01);
    end
end
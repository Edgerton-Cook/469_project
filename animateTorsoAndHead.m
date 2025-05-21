% function animateTorsoAndHead(t, theta_t, theta_h, lt, lh)
% %ANIMATETORSOANDHEAD Animates a torso with attached head
% %   Displays an animation of a rod (torso) with another rod (head) attached
% %   
% %   Input:
% %       t       - time vector
% %       theta_t - torso angle vector in degrees
% %       theta_h - head angle vector in degrees
% %       lt      - length of the torso rod
% %       lh      - length of the head rod

%     % Create figure
%     figure('Position', [100, 100, 800, 600]);
%     h_ax = axes('XLim', [-1.5*lt, 1.5*lt], 'YLim', [-0.2, 1.5*lt]);
%     axis equal;
%     grid on;
%     hold on;

%     % Wall on the left side with hatching
%     wall_x = -1*lt;
%     wall = line([wall_x, wall_x], [0, lt], 'LineWidth', 1, 'Color', 'k');
%     % Add hatching
%     for i = 0:0.05:lt
%         line([wall_x, wall_x-0.05], [i, i+0.05], 'LineWidth', 0.5, 'Color', 'k');
%     end
    
%     % Fixed pivot at bottom
%     base = line([-0.2, 0.2], [0, 0], 'LineWidth', 1, 'Color', 'k');
%     pivot = plot(0, 0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

%     % Initialize graphics objects
%     torso = line([0, 0], [0, lt], 'LineWidth', 1.5, 'Color', 'k');
%     head = line([0, 0], [lt, lt+lh], 'LineWidth', 1.2, 'Color', 'k');
    
%     % Spring connection points on wall
%     spring_wall_y = 0.8*lt;
%     spring_damper_conn = plot(wall_x, spring_wall_y, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
    
%     % Spring and damper
%     spring = line([0, 0], [0, 0], 'LineWidth', 1, 'Color', 'k');
%     damper = line([0, 0], [0, 0], 'LineWidth', 1, 'Color', 'k', 'LineStyle', ':');

%     % Add labels
%     title('Torso and Head Motion Simulation', 'FontWeight', 'normal');
%     xlabel('Position (m)');
%     ylabel('Height (m)');
    
%     % Legend
%     %legend('Wall', 'Base', 'Pivot', 'Torso', 'Head', 'Connection', 'Spring', 'Damper', ...
%     %       'Location', 'northeast');
    
%     % Animation loop
%     for i = 1:10:length(t)
%         % Get current angles
%         torso_theta = theta_t(i) * pi/180;  % Convert to radians
%         head_theta = theta_h(i) * pi/180;   % Convert to radians
        
%         % Calculate torso endpoint
%         x_torso = lt * sin(torso_theta);
%         y_torso = lt * cos(torso_theta);
        
%         % Calculate head endpoint relative to torso endpoint
%         % The head angle is relative to the torso
%         x_head = x_torso + lh * sin(torso_theta + head_theta);
%         y_head = y_torso + lh * cos(torso_theta + head_theta);
        
%         % Update torso position
%         set(torso, 'XData', [0, x_torso], 'YData', [0, y_torso]);
        
%         % Update head position
%         set(head, 'XData', [x_torso, x_head], 'YData', [y_torso, y_head]);
        
%         % Calculate spring attachment point (1/3 up the torso)
%         spring_torso_x = x_torso/3;
%         spring_torso_y = y_torso/3;
        
%         % Update spring position
%         set(spring, 'XData', [wall_x, spring_torso_x], 'YData', [spring_wall_y, spring_torso_y]);
        
%         % Update damper position (slightly lower than spring)
%         set(damper, 'XData', [wall_x, spring_torso_x], 'YData', [spring_wall_y-0.1, spring_torso_y-0.05]);
        
%         % Draw angle indicators
%         delete(findobj(h_ax, 'Tag', 'angle'));
        
%         % Torso angle arc
%         r_torso = 0.15;
%         th_torso = linspace(0, torso_theta, 20);
%         x_arc_torso = r_torso * sin(th_torso);
%         y_arc_torso = r_torso * cos(th_torso);
%         plot(x_arc_torso, y_arc_torso, 'k-', 'LineWidth', 1, 'Tag', 'angle');
        
%         % Head angle arc (relative to torso)
%         r_head = 0.1;
%         th_head = linspace(0, head_theta, 20);
%         x_arc_head = x_torso + r_head * sin(torso_theta + th_head);
%         y_arc_head = y_torso + r_head * cos(torso_theta + th_head);
%         plot([x_torso, x_arc_head], [y_torso, y_arc_head], 'k-', 'LineWidth', 1, 'Tag', 'angle');
        
%         % Add circles at each joint
%         % Head-torso joint
%         head_joint = plot(x_torso, y_torso, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'Tag', 'joint');
        
%         % Add text for current time and angles
%         delete(findobj(h_ax, 'Tag', 'info_text'));
%         text(-1.4*lt, 1.4*lt, sprintf('Time: %.2f s', t(i)), 'FontSize', 10, 'Tag', 'info_text');
%         text(-1.4*lt, 1.3*lt, sprintf('Torso: %.2f deg', theta_t(i)), 'FontSize', 10, 'Tag', 'info_text');
%         text(-1.4*lt, 1.2*lt, sprintf('Head: %.2f deg', theta_h(i)), 'FontSize', 10, 'Tag', 'info_text');
        
%         drawnow;
%         pause(0.01);
%     end

%     % Clean up labels if needed
%     delete(findobj(h_ax, 'Tag', 'info_text'));
%     delete(findobj(h_ax, 'Tag', 'angle'));
%     delete(findobj(h_ax, 'Tag', 'joint'));
% end

function animateTorsoAndHead(t, theta_t, theta_h, lt, lh)
%ANIMATETORSOANDHEAD Animates a torso with attached head
%   Displays an animation of a rod (torso) with another rod (head) attached
%   
%   Input:
%       t       - time vector
%       theta_t - torso angle vector in degrees
%       theta_h - head angle vector in degrees
%       lt      - length of the torso rod
%       lh      - length of the head rod

    % Create figure
    figure('Position', [100, 100, 800, 600]);
    h_ax = axes('XLim', [-1.5*lt, 1.5*lt], 'YLim', [-0.2, 1.5*lt]);
    axis equal;
    grid on;
    hold on;

    % Wall on the left side with hatching
    wall_x = -1*lt;
    wall = line([wall_x, wall_x], [0, lt], 'LineWidth', 1, 'Color', 'k');
    % Add hatching
    for i = 0:0.05:lt
        line([wall_x, wall_x-0.05], [i, i+0.05], 'LineWidth', 0.5, 'Color', 'k');
    end
    
    % Fixed pivot at bottom
    base = line([-0.2, 0.2], [0, 0], 'LineWidth', 1, 'Color', 'k');
    pivot = plot(0, 0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

    % Initialize graphics objects
    torso = line([0, 0], [0, lt], 'LineWidth', 1.5, 'Color', 'k');
    head = line([0, 0], [lt, lt+lh], 'LineWidth', 1.2, 'Color', 'k');
    
    % Spring connection points on wall
    spring_wall_y = 0.8*lt;
    spring_damper_conn = plot(wall_x, spring_wall_y, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
    
    % Spring and damper
    spring = line([0, 0], [0, 0], 'LineWidth', 1, 'Color', 'k');
    damper = line([0, 0], [0, 0], 'LineWidth', 1, 'Color', 'k');

    % Add labels
    title('Torso and Head Motion Simulation', 'FontWeight', 'normal');
    xlabel('Position (m)');
    ylabel('Height (m)');
    
    % Animation loop
    for i = 1:10:length(t)
        % Get current angles
        torso_theta = theta_t(i) * pi/180;  % Convert to radians
        head_theta = theta_h(i) * pi/180;   % Convert to radians
        
        % Calculate torso endpoint
        x_torso = lt * sin(torso_theta);
        y_torso = lt * cos(torso_theta);
        
        % Calculate head endpoint relative to torso endpoint
        % The head angle is relative to the torso
        x_head = x_torso + lh * sin(torso_theta + head_theta);
        y_head = y_torso + lh * cos(torso_theta + head_theta);
        
        % Update torso position
        set(torso, 'XData', [0, x_torso], 'YData', [0, y_torso]);
        
        % Update head position
        set(head, 'XData', [x_torso, x_head], 'YData', [y_torso, y_head]);
        
        % Calculate spring attachment point (1/3 up the torso)
        spring_torso_x = x_torso/3;
        spring_torso_y = y_torso/3;
        
        % Update spring position
        set(spring, 'XData', [wall_x, spring_torso_x], 'YData', [spring_wall_y, spring_torso_y]);
        
        % Update damper position (slightly lower than spring)
        set(damper, 'XData', [wall_x, spring_torso_x], 'YData', [spring_wall_y-0.1, spring_torso_y-0.05]);
        
        % Add head-torso joint
        head_joint = plot(x_torso, y_torso, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'w');
        
        % Add text for current time and angles
        delete(findobj(h_ax, 'Tag', 'info_text'));
        text(-1.4*lt, 1.4*lt, sprintf('Time: %.2f s', t(i)), 'FontSize', 10, 'Tag', 'info_text');
        text(-1.4*lt, 1.3*lt, sprintf('Torso: %.2f deg', theta_t(i)), 'FontSize', 10, 'Tag', 'info_text');
        text(-1.4*lt, 1.2*lt, sprintf('Head: %.2f deg', theta_h(i)), 'FontSize', 10, 'Tag', 'info_text');
        
        drawnow;
        pause(0.01);
        
        % Delete temporary elements
        if exist('head_joint', 'var') && ishandle(head_joint)
            delete(head_joint);
        end
    end
    
    % Clean up text at the end
    delete(findobj(h_ax, 'Tag', 'info_text'));
end
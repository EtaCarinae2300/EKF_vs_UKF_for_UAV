% Name - Surya Shashank Sekhar, Hari
% Date - Jan 14, 2026
% Department - Aerospace Engineering Sciences
% Course - ASEN 5128 Small UAS

%% Plot Simulation

function PlotSimulation(time, aircraft_state_array, control_input_array, col)
    
    % 1. Data Extraction (Vectorized)
    % -------------------------------------------------------------------
    % Extract state blocks
    pos_inertial  = aircraft_state_array(:, 1:3);   % [x, y, z]
    euler_angles  = aircraft_state_array(:, 4:6);   % [phi, theta, psi]
    vel_body      = aircraft_state_array(:, 7:9);   % [u, v, w]
    rates_body    = aircraft_state_array(:, 10:12); % [p, q, r]
    
    % Derived quantities for plotting
    h = -pos_inertial(:, 3);              % Height is -z
    euler_deg = rad2deg(euler_angles);    % Convert to degrees
    rates_deg = rad2deg(rates_body);      % Convert to deg/s
    
    % Extract controls
    delta_e = rad2deg(control_input_array(:, 1));
    delta_a = rad2deg(control_input_array(:, 2));
    delta_r = rad2deg(control_input_array(:, 3));
    delta_t = control_input_array(:, 4);  % Throttle is unitless
    
% Figure 1: States
    figure(1); 
    
    % --- Row 1: Inertial Position (x, y, z) ---
    subplot(4, 3, 1); hold on; plot(time, pos_inertial(:,1), col); ylabel('North [m]'); grid on;
    subplot(4, 3, 2); hold on; plot(time, pos_inertial(:,2), col); ylabel('East [m]');  grid on; title('Inertial Position');
    subplot(4, 3, 3); hold on; plot(time, pos_inertial(:,3), col); ylabel('Down [m]');  grid on;
    
    % --- Row 2: Euler Angles (phi, theta, psi) ---
    subplot(4, 3, 4); hold on; plot(time, euler_deg(:,1), col); ylabel('\phi [deg]');   grid on;
    subplot(4, 3, 5); hold on; plot(time, euler_deg(:,2), col); ylabel('\theta [deg]'); grid on; title('Euler Angles');
    subplot(4, 3, 6); hold on; plot(time, euler_deg(:,3), col); ylabel('\psi [deg]');   grid on;
    
    % --- Row 3: Body Velocities (u, v, w) ---
    subplot(4, 3, 7); hold on; plot(time, vel_body(:,1), col); ylabel('u [m/s]'); grid on;
    subplot(4, 3, 8); hold on; plot(time, vel_body(:,2), col); ylabel('v [m/s]'); grid on; title('Body Velocities');
    subplot(4, 3, 9); hold on; plot(time, vel_body(:,3), col); ylabel('w [m/s]'); grid on;
    
    % --- Row 4: Angular Rates (p, q, r) ---
    subplot(4, 3, 10); hold on; plot(time, rates_deg(:,1), col); ylabel('p [deg/s]'); grid on;
    subplot(4, 3, 11); hold on; plot(time, rates_deg(:,2), col); ylabel('q [deg/s]'); grid on; title('Angular Rates');
    subplot(4, 3, 12); hold on; plot(time, rates_deg(:,3), col); ylabel('r [deg/s]'); grid on;
    
    % Figure 2: Controls
    figure(2); 
    subplot(4, 1, 1); hold on; plot(time, delta_e, col); ylabel('\delta_e [deg]'); grid on; title('Control Inputs');
    subplot(4, 1, 2); hold on; plot(time, delta_a, col); ylabel('\delta_a [deg]'); grid on;
    subplot(4, 1, 3); hold on; plot(time, delta_r, col); ylabel('\delta_r [deg]'); grid on;
    subplot(4, 1, 4); hold on; plot(time, delta_t, col); ylabel('\delta_t [0-1]'); grid on; xlabel('Time [s]');
    
    % 7. Figure 6: 3D Path (Positive Height Upward)
    % -------------------------------------------------------------------
    figure(3); hold on;
    % Plot trajectory
    plot3(pos_inertial(:,2), pos_inertial(:,1), h, col, 'LineWidth', 1.5);
    
    % Start Marker (Green Circle)
    plot3(pos_inertial(1,2), pos_inertial(1,1), h(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    
    % Finish Marker (Red Square)
    plot3(pos_inertial(end,2), pos_inertial(end,1), h(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    
    grid on; axis equal;
    xlabel('East [m]'); ylabel('North [m]'); zlabel('Height [m]');
    title('3D Path');
    legend('Trajectory', 'Start', 'Finish');

end

% Name - Surya Shashank Sekhar, Hari
% Date - Jan 14, 2026
% Department - Aerospace Engineering Sciences
% Course - ASEN 5128 Small UAS

%% Aircraft Equations of Motion

function [xdot] = AircraftEOM(time,aircraft_state,aircraft_surfaces,wind_inertial,aircraft_parameters)
    
    % Aircraft State variables : 12x1 vector [x y z euler_angles(1) euler_angles(2) psi u v w p q r]
    position_inertial = aircraft_state(1:3); %[x_E, y_E, z_E]
    euler_angles = aircraft_state(4:6); % [euler_angles(1); euler_angles(2); psi]
    velocity_body = aircraft_state(7:9); % [u^E; v^E; w^E]
    rates_body    = aircraft_state(10:12); % [p; q; r]
    
    % Initializing other variables
    m = aircraft_parameters.m;
    [density, ~, ~, ~, ~] = stdatmo(-position_inertial(3));
    Transform_matrix = [1 sin(euler_angles(1))*tan(euler_angles(2)) cos(euler_angles(1))*tan(euler_angles(2)); 0 cos(euler_angles(1)) -sin(euler_angles(1)); 0 sin(euler_angles(1))*sec(euler_angles(2)) cos(euler_angles(1))*sec(euler_angles(2))];

    [aircraft_forces, aircraft_moments] = AircraftForcesAndMoments(aircraft_state, aircraft_surfaces, wind_inertial, density, aircraft_parameters);
    
    xdot = zeros(12, 1);
    xdot(1:3) = TransformFromBodyToInertial(velocity_body, euler_angles); % Position derivatives
    xdot(4:6) = Transform_matrix*rates_body;    % Euler angle derivatives
    xdot(7:9) = -cross(rates_body,velocity_body) + aircraft_forces / m; % Velocity derivatives
    xdot(10:12) = inv(aircraft_parameters.inertia_matrix) * (-cross(rates_body, (aircraft_parameters.inertia_matrix*rates_body)) + aircraft_moments); % Moment derivatives
end

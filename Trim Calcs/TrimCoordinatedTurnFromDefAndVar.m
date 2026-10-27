function [aircraft_state, control_surfaces] = TrimCoordinatedTurnFromDefAndVar(trim_definition, trim_variables)
    % Definitions: [Va, gamma, h, R] 
    Va = trim_definition(1) ; gamma = trim_definition(2); 
    h = trim_definition(3)  ; R = trim_definition(4);
    
    % Variables: [alpha, de, dt, phi, beta, da, dr]
    alpha = trim_variables(1) ; de = trim_variables(2); 
    dt = trim_variables(3)    ; phi = trim_variables(4); 
    beta = trim_variables(5)  ; da = trim_variables(6); 
    dr = trim_variables(7);

    % Position
    pos_inertial = [0; 0; -h];
    
    % Euler Angles
    euler_angles = [phi; alpha + gamma; 0]; % theta = alpha + gamma;

    % Velocity 
    vel_body = WindAnglesToAirRelativeVelocityVector([Va; beta; alpha]);

    % Rotation Rates 
    psi_dot = (Va * cos(gamma)) / R;
    rates_body = TransformFromInertialToBody([0; 0; psi_dot], euler_angles);

    % 4. Assemble Outputs
    aircraft_state = [pos_inertial; euler_angles; vel_body; rates_body];
    control_surfaces = [de; da; dr; dt];
end
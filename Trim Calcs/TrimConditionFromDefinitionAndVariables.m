function [aircraft_state, control_surfaces] = TrimConditionFromDefinitionAndVariables(trim_definition, trim_variables)
    % Trim Definition: [Va; gamma; h]
    Va    = trim_definition(1);
    gamma = trim_definition(2);
    h     = trim_definition(3);
    
    % Trim Variables: [alpha; delta_e; delta_t]
    alpha   = trim_variables(1);
    delta_e = trim_variables(2);
    delta_t = trim_variables(3);
    
    % Inertial Position Vector [x; y; z]
    % x & y do not matter so I set it to '0'
    pos_inertial = [0; 0; -h];
    
    % Euler Angles Vector [phi; theta; psi]
    % phi is defined '0' and psi does not matter
    euler_angles = [0; (alpha + gamma); 0];
    
    % Velocity Vector [u; v; w]
    % v is set to '0', rest are calculated as defined
    wind_angles = [Va; 0; alpha]; 
    vel_body = WindAnglesToAirRelativeVelocityVector(wind_angles);
    
    % Rotation Vector [p; q; r]
    % Rates are '0'
    rates_body = [0; 0; 0];
    
    % Aircraft State Vector & Control Surfaces Vector
    aircraft_state = [pos_inertial; euler_angles; vel_body; rates_body];
    control_surfaces = [delta_e; 0; 0; delta_t];
end
function [wind_inertial, v_air_body, v_air_inertial] = CalculateInertialWind(velocity_inertial, euler_angles, angular_rates, wind_angles)
    % Inputs:
    % velocity_inertial: n x 3 matrix [u_E, v_E, w_E]
    % euler_angles: n x 3 matrix [roll, pitch, yaw] (in DEGREES)
    % angular_rates: n x 3 matrix [p, q, r] (in DEGREES/S) 
    % wind_angles: n x 3 matrix [Va, beta, alpha] (angles in DEGREES)
    
    n = size(velocity_inertial, 1);
    
    % Preallocate outputs for speed
    wind_inertial = zeros(n, 3);
    v_air_body = zeros(n, 3);
    v_air_inertial = zeros(n, 3);
    
    for i = 1:n
        % 1. Extract and convert degrees to radians for your functions
        Va = wind_angles(i, 1);
        beta_rad = deg2rad(wind_angles(i, 2));
        alpha_rad = deg2rad(wind_angles(i, 3));
        euler_rad = deg2rad(euler_angles(i, :)');
        
        % 2. Air relative velocity in body frame
        v_air_body(i, :) = WindAnglesToAirRelativeVelocityVector([Va; beta_rad; alpha_rad])';
        
        % 3. Air relative velocity in inertial frame
        v_air_inertial(i, :) = TransformFromBodyToInertial(v_air_body(i, :)', euler_rad)';
        
        % 4. Inertial wind velocity
        % W_inertial = V_aircraft_inertial - V_air_inertial
        % (Ignoring sensor tip offset per instructions)
        wind_inertial(i, :) = velocity_inertial(i, :) - v_air_inertial(i, :);
    end
end
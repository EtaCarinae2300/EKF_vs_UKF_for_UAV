function cost = CostForCoordinatedTurn(trim_variables, trim_definition, aircraft_parameters)
    % State and Control Vectors
    [aircraft_state, control_surfaces] = TrimCoordinatedTurnFromDefAndVar(trim_definition, trim_variables);
    
    % Calculate Actual State Derivatives (f_EOM)
    f_EOM = AircraftEOM(0, aircraft_state, control_surfaces, [0;0;0], aircraft_parameters); 
    
    % Calculate Side Force (Y)  
    density = stdatmo(trim_definition(3)); % Height is in definition
    [aero_forces, ~] = AeroForcesAndMoments_BodyState_WindCoeffs(aircraft_state, control_surfaces, [0;0;0], density, aircraft_parameters);
    Y = aero_forces(2);
    
    % Define Target Trim Derivatives (x_dot_theoretical)
    Va = trim_definition(1);
    gamma = trim_definition(2);
    R = trim_definition(4);
    psi_dot = (Va * cos(gamma)) / R;
    
    % Theoretical State Derivative
    x_dot_theoretical = [0; 0; psi_dot; 0; 0; 0; 0; 0; 0];
    
    % Error Calculation
    error_vector = x_dot_theoretical - f_EOM(4:12);
    
    % Cost Function (Equation 10)
    cost = norm(error_vector)^2 + Y^2;
end
function cost = CostForStraightTrim(trim_variables, trim_definition, aircraft_parameters)
    % State and Control Vectors
    [aircraft_state, control_surfaces] = TrimConditionFromDefinitionAndVariables(trim_definition, trim_variables);
    
    % Calculate Actual State Derivatives (f_EOM)
    f_EOM = AircraftEOM(0, aircraft_state, control_surfaces, [0;0;0], aircraft_parameters); % time doesn't matter & wind is neglected for now
    
    % Theoretical Trim Derivative (x_dot_star)
    x_dot_theoretical = zeros(9, 1); % Only care about states 4:12
    
    % Error Calculation
    % e_trim = x_dot_star - f_EOM
    error_vector = x_dot_theoretical - f_EOM(4:12);
    
    % Cost Function
    cost = norm(error_vector)^2; 
end
function [aircraft_state_trim, control_surface_trim, optimal_trim_vars, fval] = CalculateTrim(trim_definition, aircraft_parameters)

    % Initial Guess for Optimizer [alpha0; de0; dt0]
    x0 = [0; 0; 0.2];
    
    % Bounds - Inputs for fmincon
    max_angle = 45*pi/180;
    lower_bound = [-max_angle; -max_angle; 0];
    upper_bound = [max_angle; max_angle; 1];
    
    % Cost Function Definition to make the code easier to read
    cost_func = @(x) CostForStraightTrim(x, trim_definition, aircraft_parameters);
    
    % EOM Calculation using fmincon
    [optimal_trim_vars, fval] = fmincon(cost_func, x0, [], [], [], [], lower_bound, upper_bound);
    [aircraft_state_trim, control_surface_trim] = TrimConditionFromDefinitionAndVariables(trim_definition, optimal_trim_vars);

end
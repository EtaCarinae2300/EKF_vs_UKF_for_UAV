function [aircraft_state_trim, control_surface_trim, fval] = CalculateCoordinatedTurn(trim_definition, aircraft_parameters)
    % Initial Guess for Optimizer [alpha; de; dt; phi; beta; da; dr]
    x0 = [0; 0; 0.5; 0.1; 0; 0; 0];
    
    % Bounds - Inputs for fmincon
    max_angle = 45*pi/180;
    lower_bound = [-max_angle; -max_angle; 0; -max_angle; -max_angle; -max_angle; -max_angle];
    upper_bound = [max_angle; max_angle; 1; max_angle; max_angle; max_angle; max_angle];
    
    % Cost Function Definition
    cost_func = @(x) CostForCoordinatedTurn(x, trim_definition, aircraft_parameters);
    
    % EOM Calculation using fmincon
    [optimal_trim_vars, fval] = fmincon(cost_func, x0, [], [], [], [], lower_bound, upper_bound);
    [aircraft_state_trim, control_surface_trim] = TrimCoordinatedTurnFromDefAndVar(trim_definition, optimal_trim_vars);
end
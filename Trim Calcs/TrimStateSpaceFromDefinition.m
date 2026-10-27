function [Alon, Blon, Alat, Blat] = TrimStateSpaceFromDefinition(trim_definition, aircraft_parameters)
    % TrimStateSpaceFromDefinition - Calculates Linear State Space Matrices
    % Analytically (No Finite Difference/Epsilon)

    % 1. Unpack Inputs
    Va_star    = trim_definition(1);
    gamma_star = trim_definition(2);
    h_star     = trim_definition(3);
    P = aircraft_parameters;

    % 2. Get Trim Condition
    [rho, ~, ~, ~] = stdatmo(h_star);
    [trim_state, trim_input] = CalculateTrim(trim_definition, P);
    
    alpha_star   = trim_state(5) - gamma_star; 
    theta_star   = trim_state(5); 
    delta_e_star = trim_input(1);
    delta_t_star = trim_input(4);

    % =========================================================================
    % LONGITUDINAL MODEL (Alon, Blon)
    % =========================================================================
    
    % --- Aerodynamic Coefficients ---
    CL_star = P.CL0 + P.CLalpha * alpha_star;
    CD_star = P.CDmin + P.K * (CL_star - P.CLmin)^2;
    CX_star = -cos(alpha_star)*CD_star + sin(alpha_star)*CL_star;
    CZ_star = -sin(alpha_star)*CD_star - cos(alpha_star)*CL_star;

    % --- Stability Derivatives (w.r.t Alpha) ---
    dCD_dAlpha = (2 * P.K * (CL_star - P.CLmin)) * P.CLalpha;
    CX_alpha = -(-sin(alpha_star)*CD_star + cos(alpha_star)*dCD_dAlpha) + ...
                (cos(alpha_star)*CL_star + sin(alpha_star)*P.CLalpha);
    CZ_alpha = -(cos(alpha_star)*CD_star + sin(alpha_star)*dCD_dAlpha) - ...
                (-sin(alpha_star)*CL_star + cos(alpha_star)*P.CLalpha);

    % --- Rate Derivatives (w.r.t q) ---
    CX_q = P.CLq*sin(alpha_star); 
    CZ_q = -P.CLq*cos(alpha_star);

    % --- Dimensional Constants ---
    QS  = 0.5 * rho * Va_star^2 * P.S; 
    
    % --- Xu Calculation (Aero + Analytical Propulsion Damping) ---
    Xu_aero = (rho * P.S * Va_star / P.m) * CX_star;
    
    % ANALYTICAL PROPULSION DERIVATIVE (dThrust/dVa)
    k = P.kmotor;
    C_prop_const = rho * P.Sprop * P.Cprop;
    
    % Derivative derived from Product Rule on Thrust Equation
    dT_dV = C_prop_const * delta_t_star * ( k - 2*Va_star + 2*delta_t_star*(Va_star - k) );
    
    Xu_prop = dT_dV / P.m;
    Xu = Xu_aero + Xu_prop;

    % --- Other Dimensional Derivatives ---
    Xalpha = (QS / P.m) * CX_alpha;
    Xq     = (QS * P.c / (2 * P.m * Va_star)) * CX_q; 
    Zu     = (rho * P.S * Va_star / P.m) * CZ_star;
    Zalpha = (QS / P.m) * CZ_alpha;
    Zq     = (QS * P.c / (2 * P.m * Va_star)) * CZ_q;
    Mu     = 0; 
    Malpha = (QS * P.c / P.Iy) * P.Cmalpha;
    Mq     = (QS * P.c^2 / (2 * P.Iy * Va_star)) * P.Cmq;

    % --- Control Derivatives (Analytic) ---
    % Elevator
    CX_de = P.CLde * sin(alpha_star);
    CZ_de = -P.CLde * cos(alpha_star);
    X_de = (QS / P.m) * CX_de;
    Z_de = (QS / P.m) * CZ_de;
    M_de = (QS * P.c / P.Iy) * P.Cmde;
    
    % ANALYTICAL THROTTLE DERIVATIVE (dThrust/dDelta_t)
    % Derivative derived from Thrust Equation w.r.t delta_t
    dT_ddt = C_prop_const * (k - Va_star) * ( Va_star + 2*delta_t_star*(k - Va_star) );
    
    X_dt = dT_ddt / P.m;

    % --- Assemble Alon (5x5) ---
    Alon = zeros(5,5);
    Alon(1,1) = Xu;
    Alon(1,2) = Xalpha;
    Alon(1,3) = Xq - (Va_star * sin(alpha_star)); 
    Alon(1,4) = -P.g * cos(theta_star);
    
    Alon(2,1) = Zu / Va_star;
    Alon(2,2) = Zalpha / Va_star;
    Alon(2,3) = 1 + (Zq / Va_star); 
    Alon(2,4) = -P.g * sin(theta_star) / Va_star;
    
    Alon(3,1) = Mu;
    Alon(3,2) = Malpha;
    Alon(3,3) = Mq;
    
    Alon(4,3) = 1;
    
    Alon(5,1) = sin(theta_star);
    Alon(5,2) = -Va_star * cos(theta_star); 
    Alon(5,4) = Va_star * cos(theta_star);

    % --- Assemble Blon (5x2) ---
    Blon = zeros(5,2);
    Blon(1,1) = X_de;           Blon(1,2) = X_dt;
    Blon(2,1) = Z_de / Va_star; Blon(2,2) = 0;
    Blon(3,1) = M_de;           Blon(3,2) = 0;
    
    % =========================================================================
    % LATERAL MODEL (Alat, Blat)
    % =========================================================================
    
    % --- Dimensional Derivatives ---
    Yv      = (rho * P.S * Va_star / (2 * P.m)) * P.CYbeta; 
    Yp      = (rho * P.S * P.b * Va_star / (4 * P.m)) * P.CYp;
    Yr      = (rho * P.S * P.b * Va_star / (4 * P.m)) * P.CYr;
    Ydelta_a = (rho * P.S * Va_star^2 / (2 * P.m)) * P.CYda;
    Ydelta_r = (rho * P.S * Va_star^2 / (2 * P.m)) * P.CYdr;

    L_beta_dim = (0.5 * rho * Va_star^2 * P.S * P.b) * P.Cpbeta;
    N_beta_dim = (0.5 * rho * Va_star^2 * P.S * P.b) * P.Crbeta;
    L_p_dim = (0.5 * rho * Va_star * P.S * P.b^2 / 2) * P.Cpp;
    L_r_dim = (0.5 * rho * Va_star * P.S * P.b^2 / 2) * P.Cpr;
    N_p_dim = (0.5 * rho * Va_star * P.S * P.b^2 / 2) * P.Crp;
    N_r_dim = (0.5 * rho * Va_star * P.S * P.b^2 / 2) * P.Crr;
    
    L_da_dim = (0.5 * rho * Va_star^2 * P.S * P.b) * P.Cpda;
    L_dr_dim = (0.5 * rho * Va_star^2 * P.S * P.b) * P.Cpdr;
    N_da_dim = (0.5 * rho * Va_star^2 * P.S * P.b) * P.Crda;
    N_dr_dim = (0.5 * rho * Va_star^2 * P.S * P.b) * P.Crdr;

    % --- Assemble Alat (5x5) ---
    Alat = zeros(5,5);
    
    Alat(1,1) = Yv; 
    Alat(1,2) = (Yp / Va_star) + sin(alpha_star);
    Alat(1,3) = (Yr / Va_star) - cos(alpha_star);
    Alat(1,4) = P.g * cos(theta_star) / Va_star;
    
    Alat(2,1) = L_beta_dim; Alat(2,2) = L_p_dim; Alat(2,3) = L_r_dim;
    Alat(3,1) = N_beta_dim; Alat(3,2) = N_p_dim; Alat(3,3) = N_r_dim;
    Alat(4,2) = 1;          Alat(4,3) = tan(theta_star);
    Alat(5,3) = 1 / cos(theta_star);

    % --- Assemble Blat (5x2) ---
    Blat = zeros(5,2);
    Blat(1,1) = Ydelta_a / Va_star; Blat(1,2) = Ydelta_r / Va_star;
    Blat(2,1) = L_da_dim;           Blat(2,2) = L_dr_dim;
    Blat(3,1) = N_da_dim;           Blat(3,2) = N_dr_dim;
end
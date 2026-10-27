% Name - Surya Shashank Sekhar, Hari
% Date - Jan 14, 2026
% Department - Aerospace Engineering Sciences
% Course - ASEN 5128 Small UAS

%% Aero Force & Moment

function [aero_force, aero_moment] = AeroForcesAndMoments_BodyState_WindCoeffs(aircraft_state, aircraft_surfaces, wind_inertial, density, aircraft_parameters)

    % Aircraft State variables : 12x1 vector [x y z euler_angles(1) euler_angles(2) psi u v w p q r]
    position_inertial = aircraft_state(1:3); %[x_E, y_E, z_E]
    euler_angles = aircraft_state(4:6); % [euler_angles(1); euler_angles(2); psi]
    velocity_body = aircraft_state(7:9); % [u^E; v^E; w^E]
    rates_body    = aircraft_state(10:12); % [p; q; r]

    % Control surfaces : 4x1 vector [delta_e delta_a delta_r delta_t]
    delta_e = aircraft_surfaces(1);
    delta_a = aircraft_surfaces(2);
    delta_r = aircraft_surfaces(3);
    delta_t = aircraft_surfaces(4);
    
    % Aerocynamic Coefficients : [from ttwistor.m data]
    P = aircraft_parameters;

    % 2. Air Relative Velocity Vector
    wind_body = TransformFromInertialToBody(wind_inertial, euler_angles);
    V = velocity_body - wind_body;
    
    % 3. Wind Angles
    wind_angles = AirRelativeVelocityVectorToWindAngles(V);
    Va = wind_angles(1);
    beta = wind_angles(2);
    alpha = wind_angles(3);


    % 4. Non Dimensional Rate
    char_lengths = [P.b; P.c; P.b];
    omega_hat = (rates_body .* char_lengths) ./ (2 * Va);
    Q = 0.5 * density * Va^2;
    
    p_hat = omega_hat(1);
    q_hat = omega_hat(2);
    r_hat = omega_hat(3);


    % 5. Aerodynamic Coefficients
    CL = (P.CL0) + (P.CLalpha * alpha) + (P.CLq * q_hat) + (P.CLde * delta_e);  % Lift
    CD = (P.CDmin) + (P.K * (CL - P.CLmin)^2);  % Drag
    CT = 2 * (P.Sprop / P.S) * P.Cprop * (delta_t / Va^2) * (Va + delta_t * (P.kmotor - Va)) * (P.kmotor - Va); % Thrust acting in X
    CY = (P.CY0) + (P.CYbeta * beta) + (P.CYp * p_hat) + (P.CYr * r_hat) + (P.CYda * delta_a) + (P.CYdr * delta_r); % Coeff in Y
    Cl = (P.Cl0) + (P.Clbeta * beta) + (P.Clp * p_hat) + (P.Clr * r_hat) + (P.Clda * delta_a) + (P.Cldr * delta_r); % Coeff in L
    Cm = (P.Cm0) + (P.Cmalpha * alpha) + (P.Cmq * q_hat) + (P.Cmde * delta_e); % Coeff in M
    Cn = (P.Cn0) + (P.Cnbeta * beta) + (P.Cnp * p_hat) + (P.Cnr * r_hat) + (P.Cnda * delta_a) + (P.Cndr * delta_r); % Coeff in N
        

    % 6. Cumulative Coefficients
    C_XYZ = [ ((-CD * cos(alpha)) + (CL * sin(alpha)) + CT);
              CY;
              ((-CD * sin(alpha)) - (CL * cos(alpha))) ];
    
    moment_coeffs = [Cl; Cm; Cn];
    moment_lengths = [P.b; P.c; P.b];
    C_moment = moment_coeffs .* moment_lengths;
                    
    % 7. Dimensional Forces and Moments
    % aero_forces: 3x1 vector [X Y Z] in Body Frame
    aero_force = Q * P.S * C_XYZ;
    
    % aero_moments: 3x1 vector [L M N] in Body Frame
    aero_moment = Q * P.S * C_moment;

end

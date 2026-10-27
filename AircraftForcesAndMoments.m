% Name - Surya Shashank Sekhar, Hari
% Date - Jan 14, 2026
% Department - Aerospace Engineering Sciences
% Course - ASEN 5128 Small UAS

%% Aircraft Forces & Moments

function [aircraft_forces, aircraft_moments] = AircraftForcesAndMoments(aircraft_state, aircraft_surfaces, wind_inertial, density, aircraft_parameters)
    m = aircraft_parameters.m;
    g = aircraft_parameters.g;
    gravityforce_inertial = [0; 0; m*g];

    euler_angles = aircraft_state(4:6); % [euler_angles(1); euler_angles(2); psi]

    gravity_force = TransformFromInertialToBody(gravityforce_inertial, euler_angles); % Gravity Force in Body Coords
    
    [aero_force, aero_moment] = AeroForcesAndMoments_BodyState_WindCoeffs(aircraft_state, aircraft_surfaces, wind_inertial, density, aircraft_parameters);
    aircraft_forces = aero_force + gravity_force; 
    aircraft_moments = aero_moment;
end
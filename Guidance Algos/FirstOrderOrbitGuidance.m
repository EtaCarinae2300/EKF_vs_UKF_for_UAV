function vel = FirstOrderOrbitGuidance(t, pos, orbit_speed, orbit_radius, orbit_center, orbit_flag, orbit_gains)
    % Call OrbitGuidance to get control objectives
    % Input 'pos' is used where 'aircraft_array' is expected
    obj = OrbitGuidance(pos, orbit_speed, orbit_radius, orbit_center, orbit_flag, orbit_gains);
    
    % Extract values from the control_objectives array:
    % obj(3) = chi_c (Course Angle)
    % obj(5) = Va_c  (Airspeed Command)
    chi_c = obj(3);
    Va_c  = obj(5);
    
    % Convert polar guidance back to Cartesian velocity for the ODE solver
    % Note: dot_h_c is obj(2), which is 0 in your current OrbitGuidance
    vx = Va_c * cos(chi_c);
    vy = Va_c * sin(chi_c);
    vz = obj(2); 
    
    vel = [vx; vy; vz];
end
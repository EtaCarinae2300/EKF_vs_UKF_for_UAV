function control_objectives = OrbitGuidance(aircraft_array, orbit_speed, orbit_radius, orbit_center, orbit_flag, orbit_gains)
    
    % pos = [x; y; z]
    err_x = aircraft_array(1) - orbit_center(1);
    err_y = aircraft_array(2) - orbit_center(2);
    r = sqrt(err_x^2 + err_y^2);
    
    % Vector field components
    g_r = -orbit_gains.kr * (r - orbit_radius);
    g_theta = orbit_flag; 
    g_z = -orbit_gains.kz * (aircraft_array(3) - orbit_center(3));

    % Alpha - Speed Scaling
    alpha = 1/sqrt(g_r^2 + g_theta^2 + g_z^2);
    
    % Rotate to Cartesian
    theta = atan2(err_y, err_x);
    T_theta = [cos(theta) -sin(theta) 0;
                sin(theta) cos(theta) 0;
                0 0 1]; 
    vel = orbit_speed * alpha * T_theta * [g_r; g_theta; g_z];
    
    % Course Angle
    chi_c = atan2(vel(2), vel(1));
    
    % Feedforward Course Rate: chi_dot = V_a / R
    chi_dot_ff = orbit_flag * (orbit_speed / orbit_radius);
    
    % 3. Output objectives: [h_c; dot_h_c; chi_c; chi_dot_ff; Va_c]
    control_objectives = [-orbit_center(3); 0; chi_c; chi_dot_ff; orbit_speed];
end
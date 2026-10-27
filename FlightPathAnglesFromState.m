function [flight_angles] = FlightPathAnglesFromState(aircraft_state)
    
    u_E = aircraft_state(7); v_E = aircraft_state(8); w_E = aircraft_state(9);
    V_body = [u_E; v_E; w_E];
    V_ground = norm(V_body); % Ground Velocity of the Aircraft, Vg = |V_E|

    phi = aircraft_state(4); theta = aircraft_state(5); psi = aircraft_state(6);
    euler_angles = [phi; theta; psi];
    V_inertial = TransformFromBodyToInertial(V_body,euler_angles);
    gamma_fpa = asin(-V_inertial(3)/V_ground); % Flight Path Angle, h_dot = -z_dot = Vg*sin(gamma)
    
    chi = atan2(V_inertial(2),V_inertial(1)); % Course Angle, tan(chi) = e_dot/n_dot

    flight_angles = [V_ground; chi; gamma_fpa];

end
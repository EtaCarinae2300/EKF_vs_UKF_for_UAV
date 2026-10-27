function control_objectives = StraightLineGuidance(pos_line, dir_line, pos_ac, kpath, chi_inf, Va_c)
    % dir_line = q;
    % pos_ac = p;
    % pos_line = r;
    
    chi_q = atan2(dir_line(2), dir_line(1));
    ep = pos_ac - pos_line;
    epy = -sin(chi_q)*(ep(1)) + cos(chi_q)*(ep(2));

    chi_c_rhs = chi_inf*(2/pi)*(atan(kpath*epy));

    chi_c = chi_q - chi_c_rhs; % Commanded Course Angle
    
    ki = [0;0;1];
    n = cross(dir_line,ki)/norm(cross(dir_line,ki));
    s = ep - (dot(ep,n))*n;
    
    hc = -pos_line(3) + (sqrt(s(1)^2 + s(2)^2))*(dir_line(3)/sqrt(dir_line(1)^2 + dir_line(2)^2));

    chi_dot_c = 0; 
    h_dot_c = 0;

    control_objectives = [hc; h_dot_c; chi_c; chi_dot_c; Va_c];

end
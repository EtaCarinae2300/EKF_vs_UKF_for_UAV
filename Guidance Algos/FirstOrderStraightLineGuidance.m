function vel = FirstOrderStraightLineGuidance(t, pos, pos_line, dir_line, kpath, chi_inf, Va_c)
    
    control_objectives = StraightLineGuidance(pos_line, dir_line, pos, kpath, chi_inf, Va_c);

    vel_d = -control_objectives(2);

    % kh = 0.1;
    % z_vel = -(control_objectives(2) + kh*(control_objectives(1)*pos(3)));
    % if (z_vel^2 > control_objectives(5)^2)
    %     z_vel = -control_objectives(2);
    % end

    plane_speed = sqrt(control_objectives(5)^2 - vel_d^2);
    vel = [plane_speed*cos(control_objectives(3)); plane_speed*sin(control_objectives(3)) ; vel_d];
    
end
function xdot = KinematicGuidanceAlgorithm1a(state, control_objective, wind_inertial, params)
    % state vector: pn, pe, chi, chi_dot, h, h_dot, Va
    % control objectives: h_c, h_dot_c, chi_c, chi_dot_c, Va_c

    psi_rhs = (1/state(7))*(-wind_inertial(1)*sin(state(3)) + wind_inertial(2)*cos(state(3)));
    psi = state(3) - psi_rhs;

    pn_dot = (state(7)*cos(psi)) + wind_inertial(1);
    pe_dot = (state(7)*sin(psi)) + wind_inertial(2);

    chi_dot = state(4);
    chi_dotdot = params.bcd*(control_objective(4)-state(4)) + params.bc*(control_objective(3)-state(3));

    h_dot = state(6);
    h_dotdot = -params.ahd*control_objective(2) - params.bhd*state(6) + params.bh*(control_objective(1)-state(5));

    Va_dot = params.bva*(control_objective(5)-state(7));

    xdot = [pn_dot, pe_dot, chi_dot, chi_dotdot, h_dot, h_dotdot, Va_dot]';
end
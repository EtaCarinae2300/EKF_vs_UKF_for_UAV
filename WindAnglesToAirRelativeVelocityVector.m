% Name - Surya Shashank Sekhar, Hari
% Date - Jan 14, 2026
% Department - Aerospace Engineering Sciences
% Course - ASEN 5128 Small UAS

function [velocity_body] = WindAnglesToAirRelativeVelocityVector(wind_angles)
    velocity_body = zeros(3,1);
    Va = wind_angles(1);
    beta = wind_angles(2);
    alpha = wind_angles(3);
    velocity_body = Va * [cos(alpha)*cos(beta); sin(beta); sin(alpha)*cos(beta)];
end

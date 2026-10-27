% Name - Surya Shashank Sekhar, Hari
% Date - Jan 14, 2026
% Department - Aerospace Engineering Sciences
% Course - ASEN 5128 Small UAS

function  [wind_angles] = AirRelativeVelocityVectorToWindAngles(velocity_body)
    wind_angles = zeros(3,1);
    u = velocity_body(1);
    v = velocity_body(2);
    w = velocity_body(3);
    wind_angles(1) = sqrt(u^2 + v^2 + w^2);
    wind_angles(2) = asin(v/wind_angles(1));
    wind_angles(3) = atan2(w,u);
end
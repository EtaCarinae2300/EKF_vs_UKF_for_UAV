% Name - Surya Shashank Sekhar, Hari
% Date - Jan 14, 2026
% Department - Aerospace Engineering Sciences
% Course - ASEN 5128 Small UAS

function  R = RotationMatrix321(euler_angles)
    euler_angles(1) = euler_angles(1);
    euler_angles(2) = euler_angles(2);
    psi = euler_angles(3);
    R = [cos(euler_angles(2))*cos(psi), cos(euler_angles(2))*sin(psi), -sin(euler_angles(2));
         sin(euler_angles(1))*sin(euler_angles(2))*cos(psi)-cos(euler_angles(1))*sin(psi), sin(euler_angles(1))*sin(euler_angles(2))*sin(psi)+cos(euler_angles(1))*cos(psi), sin(euler_angles(1))*cos(euler_angles(2));
         cos(euler_angles(1))*sin(euler_angles(2))*cos(psi)+sin(euler_angles(1))*sin(psi), cos(euler_angles(1))*sin(euler_angles(2))*sin(psi)-sin(euler_angles(1))*cos(psi), cos(euler_angles(1))*cos(euler_angles(2))];
end
% Name - Surya Shashank Sekhar, Hari
% Date - Jan 14, 2026
% Department - Aerospace Engineering Sciences
% Course - ASEN 5128 Small UAS

function  vector_body = TransformFromInertialToBody(vector_inertial, euler_angles)
    R = RotationMatrix321(euler_angles);
    vector_body = R * vector_inertial; 
end
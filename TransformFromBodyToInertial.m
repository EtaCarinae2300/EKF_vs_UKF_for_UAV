% Name - Surya Shashank Sekhar, Hari
% Date - Jan 14, 2026
% Department - Aerospace Engineering Sciences
% Course - ASEN 5128 Small UAS

function  vector_inertial = TransformFromBodyToInertial(vector_body, euler_angles)
    R = RotationMatrix321(euler_angles);
    vector_inertial = R' * vector_body; 
end
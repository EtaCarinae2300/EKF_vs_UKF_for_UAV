function Magneto_sensors = MagnetometerSensor(aircraft_state, sensor_params)
    % Earth's magnetic field vector in Inertial Frame 
    % [North; East; Down]
    mag_inertial = [20544; 2842.5; 46682]; % nTesla (My Location Values)
    
    % True Euler angles
    euler = aircraft_state(4:6);
    
    % Rotate magnetic field into Body Frame
    mag_body = TransformFromInertialToBody(mag_inertial, euler);
    
    % Add sensor noise
    Magneto_sensors = mag_body + sensor_params.sig_mag * randn(3,1);
end
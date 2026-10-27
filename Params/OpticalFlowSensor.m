function OpticalFlow_sensors = OpticalFlowSensor(aircraft_state, sensor_params)
    u_true = aircraft_state(7);
    v_true = aircraft_state(8);
    h = -aircraft_state(3); 
    % 
    % % Gaussian Noise
    % u_meas = u_true + (sensor_params.sig_pixel * h) * randn(1);
    % v_meas = v_true + (sensor_params.sig_pixel * h) * randn(1);
    % 
    % Ridiculous Noise
    b = sensor_params.sig_pixel * h;
    noise_u = b * sign(rand-0.5) .* log(max(eps, 1 - 2*abs(rand-0.5)));
    noise_v = b * sign(rand-0.5) .* log(max(eps, 1 - 2*abs(rand-0.5)));

    u_meas = u_true + noise_u;
    v_meas = v_true + noise_v;

    % Simulate a 1% chance the sensor just loses lock and returns 0
    if rand < 0.01
        u_meas = 0;
        v_meas = 0;
    end
    
    OpticalFlow_sensors = [u_meas; v_meas];
end
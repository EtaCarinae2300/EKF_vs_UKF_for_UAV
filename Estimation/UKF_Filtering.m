function [aircraft_state_est, wind_inertial_est] = UKF_Filtering(time, aircraft_surfaces, aircraft_state0, OpticalFlow_sensors, IMU_sensors, Magneto_sensors, aircraft_parameters, sensor_params)

persistent xhat
persistent P

persistent press_stat
persistent press_dyn

persistent run_count
if isempty(run_count)
    run_count = 0;
end
run_count = run_count + 1;

if mod(run_count, 10) == 0
    fprintf('UKF Time Step: %d (Time: %.2f s)\n', run_count, time);
end

h_ground = sensor_params.h_ground;
density = stdatmo(h_ground);

Ts_imu = sensor_params.Ts_imu; 
g = sensor_params.g;

y_k1 = [IMU_sensors; Magneto_sensors; OpticalFlow_sensors];

n_states = 12; % Number of States
alpha = 1;
beta = 0;
kappa = 0;

lambda = ((alpha^2) * (n_states + kappa)) - n_states;

w0_m = lambda / (n_states + lambda);
w0_c = w0_m + 1 - alpha^2 + beta;

wi_m = 1 / (2*(n_states + lambda));
wi_c = wi_m;

W_m = [w0_m; repmat(wi_m, 2*n_states, 1)]; 
W_c = [w0_c; repmat(wi_c, 2*n_states, 1)]; 

if isempty(xhat)
    xhat = aircraft_state0; 
    P = eye(12) * 0.1; % Small Initial Uncertainity
end


%% Dynamics Prediction Step 

Q = diag([ ...
    1e-3 * ones(3,1); ... % Position 
    1e-6 * ones(3,1); ... % Euler angles
    1e-3 * ones(3,1); ... % Velocity
    1e-6 * ones(3,1)  ... % Angular rates
]);

% Finding Sigma Points
P = 0.5 * (P + P');
P = 0.5*(P + P') + 1e-9 * eye(n_states);
%Sk = chol(P);
try
    Sk = chol(P);
catch ME
    disp('--- UKF CRASHED @ S_k---');
    fprintf('Crashed at Time: %.2f\n', time);
    disp('Eigenvalues of P:');
    disp(eig(P));
    disp('Eigenvalues of P_minus from previous time step:');
    disp(eig(P_minus));
    rethrow(ME); % This stops the script and shows the original error
end

factor = sqrt(n_states + lambda);

Xk = zeros(n_states, 2*n_states + 1);
Xk(:,1) = xhat;

for i = 1:n_states
    Xk(:, i+1)   = xhat + factor * Sk(i, :)';
    Xk(:, n_states+i+1) = xhat - factor * Sk(i, :)'; 
end

% % Propagate each Xk through Nonlinear Dynamics (ODE45)
% f = @(t, x) AircraftEOM(t, x, aircraft_surfaces, [0;0;0], aircraft_parameters);
% 
% Xk_1 = zeros(n_states, 2*n_states + 1);
% 
% for i = 1:2*n_states + 1
%     [~, x_out] = ode45(f, [time, time + Ts_imu], Xk(:, i));
%     Xk_1(:, i) = x_out(end, :)';   % X(k+1)
% end

% % Propagate each Xk through Nonlinear Dynamics (Euler Integration)
% f = @(t, x) AircraftEOM(t, x, aircraft_surfaces, [0;0;0], aircraft_parameters);
% Xk_1 = zeros(n_states, 2*n_states + 1);
% 
% for i = 1:2*n_states + 1
%     xdot = f(time, Xk(:, i));
%     Xk_1(:, i) = Xk(:, i) + xdot * Ts_imu;   
% end

% Propagate each Xk through Nonlinear Dynamics (RK4 Integration)
f = @(t, x) AircraftEOM(t, x, aircraft_surfaces, [0;0;0], aircraft_parameters);
Xk_1 = zeros(n_states, 2*n_states + 1);

for i = 1:2*n_states + 1
    Xk_1(:, i) = RK4_Integration(f, time, Xk(:, i), Ts_imu);
end

% Recombine resultant pts to get predicted mean and covar
xhat_minus = zeros(n_states,1);
P_minus = zeros(n_states);

for i = 1:2*n_states+1
    xhat_minus = xhat_minus + W_m(i) * Xk_1(:,i);
end

for i = 1:2*n_states+1
    X_diff = Xk_1(:,i) - xhat_minus; 
    P_minus = P_minus + W_c(i) * (X_diff * X_diff'); 
end
P_minus = P_minus + Q; %% SHOULD DEFINE Q

%% Measurement Update Step

R = diag([ ...
    sensor_params.sig_accel^2 * ones(3,1); ...
    sensor_params.sig_gyro^2 * ones(3,1); ...
    sensor_params.sig_abs_press^2; ...
    sensor_params.sig_dyn_press^2; ...
    sensor_params.sig_mag^2 * ones(3,1); ...
    sensor_params.sig_pixel^2 * ones(2,1) ...
]);

% Finding New Sigma Points
P_minus = 0.5 * (P_minus + P_minus');
P_minus = 0.5*(P_minus + P_minus') + 1e-9 * eye(n_states);
%Sbar_k1 = chol(P_minus);
try
    Sbar_k1 = chol(P_minus);
catch ME
    disp('--- UKF CRASHED @ Sbar_k1 ---');
    fprintf('Crashed at Time: %.2f\n', time);
    disp('Eigenvalues of P:');
    disp(eig(P));
    disp('Eigenvalues of P_minus:');
    disp(eig(P_minus));
    rethrow(ME); % This stops the script and shows the original error
end
factor = sqrt(n_states + lambda);

X_k1 = zeros(n_states, 2*n_states + 1);
X_k1(:,1) = xhat_minus;

for i = 1:n_states
    X_k1(:, i+1)   = xhat_minus + factor * Sbar_k1(i, :)';
    X_k1(:, n_states+i+1) = xhat_minus - factor * Sbar_k1(i, :)'; 
end

% Propagate each Xk through Nonlinear Measurement Dynamics
n_meas = 13;
Gamma = zeros(n_meas, 2*n_states + 1);

for i = 1:2*n_states + 1
    Gamma(:, i) = MeasurementModel(X_k1(:, i), aircraft_surfaces, [0;0;0], aircraft_parameters, sensor_params);
end

% Recombine resultant pts to get predicted mean and covar
yhat_k1 = zeros(n_meas,1);
S_k1 = zeros(n_meas);

for i = 1:2*n_states+1
    yhat_k1 = yhat_k1 + W_m(i) * Gamma(:,i);
end

for i = 1:2*n_states+1
    Gamma_diff = Gamma(:,i) - yhat_k1; 
    S_k1 = S_k1 + W_c(i) * (Gamma_diff * Gamma_diff'); 
end
S_k1 = S_k1 + R; %% SHOULD DEFINE R

% Get state-measurement cross-covariance matrix
Cxy = zeros(n_states, n_meas);

for i = 1:2*n_states+1
    X_diff = X_k1(:,i) - xhat_minus;
    Gamma_diff = Gamma(:,i) - yhat_k1; 
    Cxy = Cxy + W_c(i)*(X_diff*Gamma_diff');
end

% Estimate Kalman gain
K_k1 = Cxy / (S_k1);

% Perform Kalman State and Covariance Update with observation y
xhat = xhat_minus + K_k1*(y_k1 - yhat_k1);
%P = P_minus - K_k1*S_k1*K_k1';
I_KH = eye(n_states) - K_k1 * (Cxy' / P_minus);
P = I_KH * P_minus * I_KH' + K_k1 * R * K_k1'; %% Joseph Form

aircraft_state_est = xhat;
wind_inertial_est = [0;0;0];

end

%% Helper Functions

function [gamma] = MeasurementModel(aircraft_state, aircraft_surfaces, wind_inertial, aircraft_parameters, sensor_params)
    
    %%% Optical Flow Model
    u_meas = aircraft_state(7);
    v_meas = aircraft_state(8);
    h = -aircraft_state(3); 

    %%% Magnetometer Model
    mag_inertial = [20544; 2842.5; 46682]; % nTesla (My Location Values)
    euler = aircraft_state(4:6);
    yaw = TransformFromInertialToBody(mag_inertial, euler);
    
    %%% IMU Model
    omega_body = aircraft_state(10:12,1);
    %euler_angles = u(4:6,1);
    roll = aircraft_state(4,1);
    pitch = aircraft_state(5,1);
    height = -aircraft_state(3,1);
    density = stdatmo(height);
    wind_body = TransformFromInertialToBody(wind_inertial, aircraft_state(4:6,1));
    air_rel_body = aircraft_state(7:9,1) - wind_body;
    wind_angles = AirRelativeVelocityVectorToWindAngles(air_rel_body);
    [fa_body, ma_body] = AeroForcesAndMoments_BodyState_WindCoeffs(aircraft_state, aircraft_surfaces, wind_inertial, density, aircraft_parameters);
    Va = wind_angles(1);
    %%% accelerometers
    y_accel = fa_body/sensor_params.m;
    %%% rate gyros
    y_gyro = omega_body;
    %%% pressure
    y_abs_pressure = density*sensor_params.g*(height-sensor_params.h_ground) + sensor_params.bias_abs_press;
    y_dyn_pressure = 0.5*density*Va*Va + sensor_params.bias_dyn_press;
    inertial_sensors = [y_accel; y_gyro; y_abs_pressure; y_dyn_pressure];
    
    gamma = [inertial_sensors; yaw; u_meas; v_meas];
end
function [aircraft_state_est, wind_inertial_est] = UKF_Filtering(time, aircraft_surfaces, aircraft_state0, OpticalFlow_sensors, IMU_sensors, Magneto_sensors, aircraft_parameters, sensor_params)

persistent xhat
persistent P

persistent run_count
if isempty(run_count)
    run_count = 0;
end
run_count = run_count + 1;

if mod(run_count, 100) == 0
    fprintf('UKF Time Step: %d (Time: %.2f s)\n', run_count, time);
end

h_ground = sensor_params.h_ground;
Ts_imu   = sensor_params.Ts_imu;
g        = sensor_params.g;

y_k1 = [IMU_sensors; Magneto_sensors; OpticalFlow_sensors];

n_states = 12;
alpha    = 1e-3;   % Smaller alpha improves numerical stability
beta     = 2;
kappa    = 0;

lambda = (alpha^2) * (n_states + kappa) - n_states;

w0_m = lambda / (n_states + lambda);
w0_c = w0_m + 1 - alpha^2 + beta;
wi_m = 1 / (2*(n_states + lambda));
wi_c = wi_m;

W_m = [w0_m; repmat(wi_m, 2*n_states, 1)];
W_c = [w0_c; repmat(wi_c, 2*n_states, 1)];

%%% Initialize on first call
if isempty(xhat)
    xhat = aircraft_state0;
    P    = diag([ ...
        1*ones(3,1);   ... % Position
        0.1*ones(3,1); ... % Euler angles
        1*ones(3,1);   ... % Velocity
        0.1*ones(3,1)  ... % Angular rates
    ]);
end

%% =========================================================
%  PROCESS NOISE
% ==========================================================
Q = diag([ ...
    1e-2 * ones(3,1); ... % Position
    1e-4 * ones(3,1); ... % Euler angles
    1e-2 * ones(3,1); ... % Velocity
    1e-4 * ones(3,1)  ... % Angular rates
]);

%% =========================================================
%  MEASUREMENT NOISE
% ==========================================================
R = diag([ ...
    sensor_params.sig_accel^2      * ones(3,1); ...
    sensor_params.sig_gyro^2       * ones(3,1); ...
    sensor_params.sig_abs_press^2; ...
    sensor_params.sig_dyn_press^2; ...
    sensor_params.sig_mag^2        * ones(3,1); ...
    sensor_params.sig_pixel^2      * ones(2,1)  ...
]);

%% =========================================================
%  PREDICTION STEP
% ==========================================================

%%% Force P symmetric and positive definite before sigma points
P    = 0.5*(P + P') + 1e-8 * eye(n_states);
P    = nearestSPD(P);

try
    Sk = chol(P);
catch ME
    warning('UKF chol(P) failed at step %d — resetting P. Error: %s', run_count, ME.message);
    P  = diag([1*ones(3,1); 0.1*ones(3,1); 1*ones(3,1); 0.1*ones(3,1)]);
    Sk = chol(P);
end

factor = sqrt(n_states + lambda);

Xk       = zeros(n_states, 2*n_states + 1);
Xk(:, 1) = xhat;
for i = 1:n_states
    Xk(:, i+1)          = xhat + factor * Sk(i,:)';
    Xk(:, n_states+i+1) = xhat - factor * Sk(i,:)';
end

%%% Euler integration propagation
f    = @(t, x) AircraftEOM(t, x, aircraft_surfaces, [0;0;0], aircraft_parameters);
Xk_1 = zeros(n_states, 2*n_states + 1);

for i = 1:2*n_states + 1
    xdot       = f(time, Xk(:, i));
    Xk_1(:, i) = Xk(:, i) + xdot * Ts_imu;
end

%%% Predicted mean and covariance
xhat_minus = zeros(n_states, 1);
for i = 1:2*n_states + 1
    xhat_minus = xhat_minus + W_m(i) * Xk_1(:, i);
end

P_minus = zeros(n_states);
for i = 1:2*n_states + 1
    d       = Xk_1(:,i) - xhat_minus;
    P_minus = P_minus + W_c(i) * (d * d');
end
P_minus = P_minus + Q;
P_minus = 0.5*(P_minus + P_minus') + 1e-8 * eye(n_states);
P_minus = nearestSPD(P_minus);

%% =========================================================
%  MEASUREMENT UPDATE STEP
% ==========================================================

try
    Sbar_k1 = chol(P_minus);
catch ME
    warning('UKF chol(P_minus) failed at step %d — resetting. Error: %s', run_count, ME.message);
    P_minus  = Q;
    Sbar_k1  = chol(P_minus);
end

X_k1       = zeros(n_states, 2*n_states + 1);
X_k1(:, 1) = xhat_minus;
for i = 1:n_states
    X_k1(:, i+1)          = xhat_minus + factor * Sbar_k1(i,:)';
    X_k1(:, n_states+i+1) = xhat_minus - factor * Sbar_k1(i,:)';
end

%%% Push sigma points through measurement model
n_meas = 13;
Gamma  = zeros(n_meas, 2*n_states + 1);
for i = 1:2*n_states + 1
    Gamma(:, i) = MeasurementModel(X_k1(:, i), aircraft_surfaces, [0;0;0], aircraft_parameters, sensor_params);
end

%%% Check for NaN/Inf in Gamma and replace with predicted mean
for i = 1:2*n_states + 1
    if any(~isfinite(Gamma(:,i)))
        warning('UKF: Non-finite measurement sigma point %d at step %d — zeroing.', i, run_count);
        Gamma(:,i) = zeros(n_meas, 1);
    end
end

%%% Predicted measurement mean and innovation covariance
yhat_k1 = zeros(n_meas, 1);
for i = 1:2*n_states + 1
    yhat_k1 = yhat_k1 + W_m(i) * Gamma(:, i);
end

S_k1 = zeros(n_meas);
for i = 1:2*n_states + 1
    dg   = Gamma(:,i) - yhat_k1;
    S_k1 = S_k1 + W_c(i) * (dg * dg');
end
S_k1 = S_k1 + R;
S_k1 = 0.5*(S_k1 + S_k1');  % Symmetrize

%%% Cross-covariance
Cxy = zeros(n_states, n_meas);
for i = 1:2*n_states + 1
    dx  = X_k1(:,i)  - xhat_minus;
    dg  = Gamma(:,i) - yhat_k1;
    Cxy = Cxy + W_c(i) * (dx * dg');
end

%%% Kalman gain
K_k1 = Cxy / S_k1;

%%% Check for bad Kalman gain
if any(~isfinite(K_k1(:)))
    warning('UKF: Non-finite Kalman gain at step %d — skipping update.', run_count);
    xhat = xhat_minus;
    P    = P_minus;
else
    %%% State update
    innovation = y_k1 - yhat_k1;
    xhat       = xhat_minus + K_k1 * innovation;

    %%% Joseph stabilized covariance update (numerically robust)
    H_approx = Cxy' / P_minus;          % Approximate H = Cxy' * P_minus^-1
    IKH      = eye(n_states) - K_k1 * H_approx;
    P        = IKH * P_minus * IKH' + K_k1 * R * K_k1';
    P        = 0.5*(P + P') + 1e-8 * eye(n_states);
end

aircraft_state_est = xhat;
wind_inertial_est  = [0; 0; 0];

end

%% =========================================================
%  MEASUREMENT MODEL
% ==========================================================
function gamma = MeasurementModel(aircraft_state, aircraft_surfaces, wind_inertial, aircraft_parameters, sensor_params)

    %%% --- Clamp altitude to valid stdatmo range ---
    h_raw = -aircraft_state(3);
    h     = max(0, min(h_raw, 84999));   % stdatmo valid: 0–84999 m

    density = stdatmo(h);
    if ~isfinite(density) || density <= 0
        density = stdatmo(sensor_params.h_ground);  % fallback to ground density
    end

    %%% Euler angles
    roll  = aircraft_state(4);
    pitch = aircraft_state(5);

    %%% Angular rates
    omega_body = aircraft_state(10:12);

    %%% Wind and airspeed
    wind_body    = TransformFromInertialToBody(wind_inertial, aircraft_state(4:6));
    air_rel_body = aircraft_state(7:9) - wind_body;
    wind_angles  = AirRelativeVelocityVectorToWindAngles(air_rel_body);
    Va           = wind_angles(1);

    %%% Aero forces
    [fa_body, ~] = AeroForcesAndMoments_BodyState_WindCoeffs( ...
        aircraft_state, aircraft_surfaces, wind_inertial, density, aircraft_parameters);

    %%% IMU model
    y_accel        = fa_body / sensor_params.m;
    y_gyro         = omega_body;
    y_abs_pressure = density * sensor_params.g * (h - sensor_params.h_ground) + sensor_params.bias_abs_press;
    y_dyn_pressure = 0.5 * density * Va^2 + sensor_params.bias_dyn_press;

    %%% Magnetometer model
    mag_inertial = [20544; 2842.5; 46682];   % nTesla
    yaw_vec      = TransformFromInertialToBody(mag_inertial, aircraft_state(4:6));

    %%% Optical flow model
    u_meas = aircraft_state(7);
    v_meas = aircraft_state(8);

    gamma = [y_accel; y_gyro; y_abs_pressure; y_dyn_pressure; yaw_vec; u_meas; v_meas];
end

%% =========================================================
%  nearestSPD: Force matrix to nearest symmetric positive definite
% ==========================================================
function Ahat = nearestSPD(A)
    % Based on Higham (1988)
    B    = 0.5*(A + A');
    [~, p] = chol(B);
    if p == 0
        Ahat = B;
        return;
    end
    [V, D] = eig(B);
    d      = diag(D);
    d      = max(d, 1e-10);   % Floor eigenvalues
    Ahat   = V * diag(d) * V';
    Ahat   = 0.5*(Ahat + Ahat');
end

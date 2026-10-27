% =========================================================
% EKF_ConsistencyTest.m
%
% Truth Model Consistency Test for EKF_Filtering
% Runs a Monte Carlo simulation and evaluates:
%   - NEES  (Normalized Estimation Error Squared)
%   - NIS   (Normalized Innovation Squared)
%
% Both statistics are compared against 95% chi-squared
% confidence bounds to judge filter consistency.
%
% Authors:  Surya Shashank Sekhar, Hari
% Course:   ASEN 5128
% =========================================================

close all; clear all; clc;

%% -------------------------------------------------------
%  User Parameters  (tune these as needed)
% --------------------------------------------------------
N_MC       = 50;          % Number of Monte Carlo runs
Tfinal     = 200;         % Simulation duration [sec]
alpha_CI   = 0.05;        % Significance level for chi-squared bounds (95% CI)

%% -------------------------------------------------------
%  Aircraft / sensor setup  (mirrors main Project.m)
% --------------------------------------------------------
ttwistor
sensor_params = SensorParametersTtwistor(aircraft_parameters);
sensor_params.sig_pixel = 0.01;
sensor_params.sig_mag   = 500;

V_trim     = 18;
h_trim     = 1805;
gamma_trim = 0;
trim_def   = [V_trim; gamma_trim; h_trim];

[aircraft_state_trim, control_input_trim, trim_variables, ~] = ...
    CalculateTrim(trim_def, aircraft_parameters);

wind_inertial  = [0; 10; 0];
Ts             = sensor_params.Ts_imu;          % 0.1 s
n_steps        = round(Tfinal / Ts);

% State indices for NEES (use all 12 states)
n_x = 12;

% Measurement dimension for NIS (optical flow u,v + baro h = 3)
n_z = 3;

%% -------------------------------------------------------
%  Pre-allocate accumulators
% --------------------------------------------------------
NEES_sum = zeros(1, n_steps);   % sum over MC runs at each step
NIS_sum  = zeros(1, n_steps);   % sum over MC runs at each step
time_vec = Ts * (0 : n_steps-1);

%% -------------------------------------------------------
%  Monte Carlo Loop
% --------------------------------------------------------
fprintf('\nRunning %d Monte Carlo trials ...\n', N_MC);

for mc = 1 : N_MC

    fprintf('  Trial %d / %d\n', mc, N_MC);

    % --- Clear EKF persistent state between runs ---
    clear EKF_Filtering

    % --- Initial true state (add small random perturbation) ---
    aircraft_state0          = aircraft_state_trim;
    aircraft_state0(3,1)     = -1655;
    aircraft_state0(4,1)     = 0;

    wind_body = TransformFromInertialToBody(wind_inertial, aircraft_state0(4:6,1));
    aircraft_state0(7:9,1)   = aircraft_state0(7:9,1) + wind_body;

    % Small random IC perturbation to spread the Monte Carlo ensemble
    ic_noise          = zeros(12,1);
    ic_noise(1:3)     = 5  * randn(3,1);         % position [m]
    ic_noise(4:6)     = deg2rad(2) * randn(3,1); % attitude [rad]
    ic_noise(7:9)     = 0.5 * randn(3,1);        % velocity [m/s]
    ic_noise(10:12)   = deg2rad(1) * randn(3,1); % rates [rad/s]
    aircraft_state0   = aircraft_state0 + ic_noise;

    control_input0 = control_input_trim;

    % --- Simulate and filter ---
    x_true = aircraft_state0;

    for i = 1 : n_steps

        TSPAN = Ts * [i-1,  i];

        % ---- Sensor measurements (true state + noise) ----
        inertial_sensors = InertialSensors(x_true, control_input0, ...
                               wind_inertial, aircraft_parameters, sensor_params);
        magneto_sensors  = MagnetometerSensor(x_true, sensor_params);
        optical_sensors  = OpticalFlowSensor(x_true, sensor_params);

        % ---- EKF update ----
        [x_est, ~, P_approx] = EKF_Filtering(time_vec(i), optical_sensors, ...
                         inertial_sensors, magneto_sensors, sensor_params);

        % ---- Retrieve covariance from EKF internals ----
        %
        % NOTE: EKF_Filtering does not currently return P.
        % You have two clean options:
        %
        %   Option A (recommended): modify EKF_Filtering to return
        %     [x_est, wind_est, P_full] where
        %       P_full = blkdiag(P_ang, P_trans)
        %     This gives exact filter covariance.
        %
        %   Option B (used here as a placeholder): build an
        %     approximate P from the sensor noise parameters.
        %     This lets the script run without touching EKF_Filtering,
        %     but NEES will only be approximate until Option A is added.
        %
        % ---- Option B placeholder covariance ----
        density = stdatmo(sensor_params.h_ground);
        g       = sensor_params.g;

        h_baro  = sensor_params.h_ground + ...
                  inertial_sensors(7) / (density * g);       % rough altitude

        sig_pos_xy  = sensor_params.sig_pixel * abs(h_baro); % [m]
        sig_pos_z   = sensor_params.sig_abs_press / (density * g);
        sig_att     = deg2rad(2);                            % 2 deg 1-sigma
        sig_vel     = 0.5;                                   % 0.5 m/s 1-sigma
        sig_rate    = deg2rad(1);                            % 1 deg/s 1-sigma
        % 
        % P_diag = [sig_pos_xy^2; sig_pos_xy^2; sig_pos_z^2;  % positions
        %           sig_att^2;   sig_att^2;   sig_att^2;       % Euler angles
        %           sig_vel^2;   sig_vel^2;   sig_vel^2;       % velocities
        %           sig_rate^2;  sig_rate^2;  sig_rate^2];     % angular rates
        % P_approx = diag(P_diag);

        % ---- NEES ----
        err              = x_true - x_est;
        % Wrap attitude angles to [-pi, pi]
        err(4:6)         = wrapToPi(err(4:6));

        NEES_i           = err' * (P_approx \ err);
        NEES_sum(i)      = NEES_sum(i) + NEES_i;

        % ---- NIS  (innovation using optical flow + baro) ----
        %
        % Predicted measurement from estimated state
        h_est_baro = h_baro;    % same baro used by filter
        u_est      = x_est(7);  % body-x velocity used as optical-flow proxy
        v_est      = x_est(8);  % body-y velocity

        z_meas     = [optical_sensors(1); optical_sensors(2); -h_baro];
        z_pred     = [u_est;              v_est;               x_est(3)];

        R_nis      = diag([sig_pos_xy^2, sig_pos_xy^2, sig_pos_z^2]);
        H_nis      = zeros(3,12);
        H_nis(1,7) = 1;  H_nis(2,8) = 1;  H_nis(3,3) = 1;
        S_nis      = H_nis * P_approx * H_nis' + R_nis;

        nu         = z_meas - z_pred;
        NIS_i      = nu' * (S_nis \ nu);
        NIS_sum(i) = NIS_sum(i) + NIS_i;

        % ---- Propagate true state ----
        [~, YOUT] = ode45(@(t,y) AircraftEOM(t, y, control_input0, ...
                              wind_inertial, aircraft_parameters), ...
                              TSPAN, x_true, []);
        x_true = YOUT(end,:)';

    end % time steps
end % Monte Carlo

fprintf('Done.\n\n');

%% -------------------------------------------------------
%  Compute averaged NEES and NIS
% --------------------------------------------------------
NEES_avg = NEES_sum / N_MC;
NIS_avg  = NIS_sum  / N_MC;

%% -------------------------------------------------------
%  Chi-squared 95% confidence bounds
%  For averaged NEES: chi2([alpha/2, 1-alpha/2], n_x*N) / N
%  For averaged NIS:  chi2([alpha/2, 1-alpha/2], n_z*N) / N
% --------------------------------------------------------
chi2_NEES_lo = chi2inv(alpha_CI/2,   n_x * N_MC) / N_MC;
chi2_NEES_hi = chi2inv(1-alpha_CI/2, n_x * N_MC) / N_MC;

chi2_NIS_lo  = chi2inv(alpha_CI/2,   n_z * N_MC) / N_MC;
chi2_NIS_hi  = chi2inv(1-alpha_CI/2, n_z * N_MC) / N_MC;

%% -------------------------------------------------------
%  Console summary
% --------------------------------------------------------
pct_NEES_inside = 100 * mean(NEES_avg >= chi2_NEES_lo & ...
                              NEES_avg <= chi2_NEES_hi);
pct_NIS_inside  = 100 * mean(NIS_avg  >= chi2_NIS_lo  & ...
                              NIS_avg  <= chi2_NIS_hi);

fprintf('====================================================\n');
fprintf('  EKF CONSISTENCY TEST RESULTS  (%d MC runs)\n', N_MC);
fprintf('====================================================\n');
fprintf('  Metric    |  Mean value  |  95%% bounds       | %% inside\n');
fprintf('  ----------+--------------+-------------------+---------\n');
fprintf('  NEES      |  %8.3f    |  [%5.2f, %5.2f]  |  %5.1f%%\n', ...
        mean(NEES_avg), chi2_NEES_lo, chi2_NEES_hi, pct_NEES_inside);
fprintf('  NIS       |  %8.3f    |  [%5.2f, %5.2f]  |  %5.1f%%\n', ...
        mean(NIS_avg),  chi2_NIS_lo,  chi2_NIS_hi,  pct_NIS_inside);
fprintf('====================================================\n');
fprintf('\n  Interpretation:\n');
fprintf('  * NEES >> upper bound  -> filter OVERCONFIDENT (P too small)\n');
fprintf('  * NEES << lower bound  -> filter UNDERCONFIDENT (P too large)\n');
fprintf('  * NIS  >> upper bound  -> innovation covariance S too small\n');
fprintf('  * ~90%%+ inside bounds  -> filter is CONSISTENT\n\n');

%% -------------------------------------------------------
%  Plots
% --------------------------------------------------------
figure('Name','EKF Consistency: NEES','NumberTitle','off');
  plot(time_vec, NEES_avg, 'b', 'LineWidth', 1.2); hold on;
  yline(chi2_NEES_lo, 'r--', 'LineWidth', 1.2, 'Label', '95% lower');
  yline(chi2_NEES_hi, 'r--', 'LineWidth', 1.2, 'Label', '95% upper');
  xlabel('Time [s]');
  ylabel(sprintf('ANEES  (n_{x}=%d, N=%d)', n_x, N_MC));
  title('Averaged NEES — EKF Consistency Test');
  legend('ANEES', '95% bounds', 'Location', 'best');
  grid on;

figure('Name','EKF Consistency: NIS','NumberTitle','off');
  plot(time_vec, NIS_avg, 'm', 'LineWidth', 1.2); hold on;
  yline(chi2_NIS_lo, 'r--', 'LineWidth', 1.2, 'Label', '95% lower');
  yline(chi2_NIS_hi, 'r--', 'LineWidth', 1.2, 'Label', '95% upper');
  xlabel('Time [s]');
  ylabel(sprintf('ANIS  (n_{z}=%d, N=%d)', n_z, N_MC));
  title('Averaged NIS — EKF Consistency Test');
  legend('ANIS', '95% bounds', 'Location', 'best');
  grid on;

figure('Name','EKF Error Statistics','NumberTitle','off');
  subplot(2,1,1);
    semilogy(time_vec, NEES_avg, 'b'); hold on;
    semilogy([time_vec(1) time_vec(end)], [chi2_NEES_lo chi2_NEES_lo], 'r--');
    semilogy([time_vec(1) time_vec(end)], [chi2_NEES_hi chi2_NEES_hi], 'r--');
    ylabel('ANEES (log)'); title('NEES — log scale'); grid on;
  subplot(2,1,2);
    semilogy(time_vec, NIS_avg, 'm'); hold on;
    semilogy([time_vec(1) time_vec(end)], [chi2_NIS_lo chi2_NIS_lo], 'r--');
    semilogy([time_vec(1) time_vec(end)], [chi2_NIS_hi chi2_NIS_hi], 'r--');
    ylabel('ANIS (log)'); xlabel('Time [s]'); title('NIS — log scale'); grid on;
function [aircraft_state_est, wind_inertial_est] = UKF_EKF_Hybrid_Filtering(time, aircraft_surfaces, aircraft_state0, OpticalFlow_sensors, IMU_sensors, Magneto_sensors, aircraft_parameters, sensor_params)
%
% magneto_sensors = [mag_x; mag_y; mag_z];
%
% OpticalFlow_sensors = [u_meas, v_meas] in Body coordinates;
%
% IMU_sensors = [y_accel; y_gyro; y_pressure; y_dyn_pressure]; 
% 

%% Variable Definition
persistent p_hat
persistent q_hat
persistent r_hat

persistent press_stat
persistent press_dyn

persistent phi_hat
persistent theta_hat
persistent psi_hat
persistent P_ang

persistent x_trans
persistent P_trans

h_ground = sensor_params.h_ground;
density = stdatmo(h_ground);

Ts_imu = sensor_params.Ts_imu; 
g = sensor_params.g;

%% Air Relative Velocity

a_Va = 0.5;%
alpha_Va = exp(-a_Va*Ts_imu);

if(isempty(press_dyn))
    press_dyn = IMU_sensors(8);
else
    press_dyn = LowPassFilter(press_dyn, IMU_sensors(8), alpha_Va);
end
Va = sqrt((2/density)*press_dyn);

%% Angular Velocities from IMU Sensor

a_omega = 1000;
alpha_omega = exp(-a_omega*Ts_imu);

if(isempty(p_hat))
    p_hat = IMU_sensors(4);
else
    p_hat = LowPassFilter(p_hat, IMU_sensors(4), alpha_omega);
end

if(isempty(q_hat))
    q_hat = IMU_sensors(5);
else
    q_hat = LowPassFilter(q_hat, IMU_sensors(5), alpha_omega);
end

if(isempty(r_hat))
    r_hat = IMU_sensors(6);
else
    r_hat = LowPassFilter(r_hat, IMU_sensors(6), alpha_omega);
end

%% Euler Angle Measurements from IMU & Magnetometer

delta = deg2rad(7.53);
if isempty(phi_hat)
    magneto_inertial = TransformFromBodyToInertial(Magneto_sensors, [0; 0; 0]);
else
    magneto_inertial = TransformFromBodyToInertial(Magneto_sensors, [phi_hat; theta_hat; 0]);
end
psi_m = -atan2(magneto_inertial(2), magneto_inertial(1));
psi_meas = delta + psi_m;

Q_scale_factor = 0.01*((pi/180)^2);
Q_ang = Q_scale_factor*eye(3);

sig_mag = 5 * pi/180;
R_ang = diag([sensor_params.sig_accel^2, sensor_params.sig_accel^2, sensor_params.sig_accel^2, sig_mag^2]);

if(isempty(phi_hat))
    phi_hat = 0;
    theta_hat = 0;
    psi_hat = psi_meas;
    P_ang = ((30*pi/180)^2)*eye(3);
else
    %%% Propagate
    [xdot_ang, A_ang] = AttitudeFilterUpdate(phi_hat, theta_hat, psi_hat, p_hat, q_hat, r_hat);
    phi_hat = phi_hat + xdot_ang(1)*Ts_imu;
    theta_hat = theta_hat + xdot_ang(2)*Ts_imu;
    psi_hat = psi_hat + xdot_ang(3)*Ts_imu;
    P_ang = P_ang + Ts_imu*(A_ang*P_ang + P_ang*A_ang' + Q_ang);
    
    %%% Measurement update
    [y_hat_ang, H_ang] = AttitudeFilterMeasurement(phi_hat, theta_hat, psi_hat, p_hat, q_hat, r_hat, Va, g);
    K_ang = (P_ang*H_ang') / (R_ang+H_ang*P_ang*H_ang');
    %P_ang_est = (eye(3)-K_ang*H_ang)*P_ang_est; % SHOULD USE (I-KH)*P*(I-KH)'+KRK'??
    P_ang = (eye(3)-K_ang*H_ang)*P_ang*(eye(3)-K_ang*H_ang)' + (K_ang*R_ang*K_ang');

    residual_ang = [IMU_sensors(1:3,1); psi_meas] - y_hat_ang;
    residual_ang(4) = wrapToPi(residual_ang(4)); 
    
    x_hat_ang = [phi_hat; theta_hat; psi_hat] + K_ang * residual_ang;

    phi_hat = x_hat_ang(1);
    theta_hat = x_hat_ang(2);
    psi_hat = x_hat_ang(3);
end



%% Translational Velocities

a_h = 1;
alpha_h = exp(-a_h*Ts_imu);

if(isempty(press_stat))
    press_stat = IMU_sensors(7);
else
    press_stat = LowPassFilter(press_stat, IMU_sensors(7), alpha_h);
end
h_hat =  h_ground + press_stat/(density*g);

%% KF on Translation

%Q_trans = eye(6); %EKF Q
Q_trans = diag([1000, 50, 1, 1, 1, 1]); % Tuned for UKF

R_xy = (sensor_params.sig_pixel * abs(h_hat))^2;
R_h = (sensor_params.sig_abs_press / (density * g))^2;
R_trans = diag([R_xy, R_xy, R_h]);
n_states_UKF = 6;

% UKF Params
lambda = ((aircraft_parameters.UKF.alpha^2) * (n_states_UKF + aircraft_parameters.UKF.kappa)) - n_states_UKF;

w0_m = lambda / (n_states_UKF + lambda);
w0_c = w0_m + 1 - aircraft_parameters.UKF.alpha^2 + aircraft_parameters.UKF.beta;

wi_m = 1 / (2*(n_states_UKF + lambda));
wi_c = wi_m;

W_m = [w0_m; repmat(wi_m, 2*n_states_UKF, 1)]; 
W_c = [w0_c; repmat(wi_c, 2*n_states_UKF, 1)];

if isempty(x_trans)
    x_trans = [0; 0; -h_ground; 0; 0; 0]; 
    P_trans = eye(6); 
else
    %%% Propagate using UKF
    P_trans = 0.5 * (P_trans + P_trans');
    Sk = chol(P_trans, 'lower');

    factor = sqrt(n_states_UKF + lambda);

    Xk = zeros(n_states_UKF, 2*n_states_UKF + 1);
    Xk(:,1) = x_trans;
    
    for i = 1:n_states_UKF
        Xk(:, i+1)   = x_trans + factor * Sk(i, :)';
        Xk(:, n_states_UKF+i+1) = x_trans - factor * Sk(i, :)'; 
    end

    % Propagate each Xk through Nonlinear Dynamics (ODE45)
    Xk_1 = zeros(n_states_UKF, 2*n_states_UKF + 1);
    accel = IMU_sensors(1:3); % Accelerometer inputs from sensor
    
    for i = 1:2*n_states_UKF + 1
        % Extract current sigma point
        sig_pt = Xk(:, i);
        
        % Calculate derivative using the kinematic helper function
        [xdot_trans, ~] = TranslationFilterUpdate(sig_pt, accel, phi_hat, theta_hat, psi_hat, p_hat, q_hat, r_hat, g);
        
        % Propagate forward using Euler integration
        Xk_1(:, i) = sig_pt + xdot_trans * Ts_imu;   
    end
    
    % Recombine resultant pts to get predicted mean and covar
    x_trans = zeros(n_states_UKF,1);
    P_trans = zeros(n_states_UKF, n_states_UKF);
    
    for i = 1:2*n_states_UKF+1
        x_trans = x_trans + W_m(i) * Xk_1(:,i);
    end
    
    for i = 1:2*n_states_UKF+1
        X_diff = Xk_1(:,i) - x_trans; 
        P_trans = P_trans + W_c(i) * (X_diff * X_diff'); 
    end
    P_trans = P_trans + Q_trans;


    %%% Measurement Update using EKF
    z_trans = [OpticalFlow_sensors(1); OpticalFlow_sensors(2); -h_hat]; % SHOULD BE NEGATIVE H? 
    [y_hat_trans, H_trans] = TranslationFilterMeasurement(x_trans);
    K_trans = (P_trans * H_trans') / (R_trans + (H_trans * P_trans * H_trans'));
    P_trans = (eye(6) - K_trans * H_trans) * P_trans;
    x_trans = x_trans + K_trans * (z_trans - y_hat_trans);
end

x_hat = x_trans(1);
y_hat = x_trans(2);
z_hat = x_trans(3);
u_hat = x_trans(4);
v_hat = x_trans(5);
w_hat = x_trans(6);

aircraft_state_est = [x_hat, y_hat, z_hat, phi_hat, theta_hat, psi_hat, u_hat, v_hat, w_hat, p_hat, q_hat, r_hat]';
wind_inertial_est = [0; 0; 0];

end

%% Helper Functions

function ynew = LowPassFilter(yold, unew, alpha)
    ynew = alpha*yold + (1-alpha)*unew;
end

function [xdot_ang, A_ang] = AttitudeFilterUpdate(phi, theta, psi, p, q, r)

xdot_ang =   [p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
          q*cos(phi) - r*sin(phi);
          q*sin(phi)*sec(theta) + r*cos(phi)*sec(theta)];
    
A_ang = [ q*cos(phi)*tan(theta) - r*sin(phi)*tan(theta),   (q*sin(phi) + r*cos(phi))/cos(theta)^2, 0;
   -q*sin(phi) - r*cos(phi), 0, 0;
   q*cos(phi)*sec(theta) - r*sin(phi)*sec(theta), q*sin(phi)*tan(theta)*sec(theta) + r*cos(phi)*tan(theta)*sec(theta), 0];

end

function [y_hat_ang, H_ang] = AttitudeFilterMeasurement(phi, theta, psi, p, q, r, Va, g)

y_accel = [q*Va*sin(theta) + g*sin(theta);
           r*Va*cos(theta) - p*Va*sin(theta) - g*cos(theta)*sin(phi);
          -q*Va*cos(theta) - g*cos(theta)*cos(phi)];

y_hat_ang = [y_accel; psi]; 

H_accel = [0, q*Va*cos(theta) + g*cos(theta), 0;
          -g*cos(phi)*cos(theta), -r*Va*sin(theta) - p*Va*cos(theta) + g*sin(phi)*sin(theta), 0;
           g*sin(phi)*cos(theta), (q*Va + g*cos(phi))*sin(theta), 0];

H_ang = [H_accel; 
     0, 0, 1]; 
end

function [xdot_trans, A_trans] = TranslationFilterUpdate(x, accel, phi, theta, psi, p, q, r, g)

    u = x(4);
    v = x(5);
    w = x(6);
    
    ax = accel(1);
    ay = accel(2);
    az = accel(3);

    R_bi = RotationMatrix321([phi, theta, psi])'; 

    p_dot = R_bi * [u; v; w];
    
    v_dot = [ax - q*w + r*v - g*sin(theta);
             ay - r*u + p*w + g*cos(theta)*sin(phi);
             az - p*v + q*u + g*cos(theta)*cos(phi)];
             
    xdot_trans = [p_dot; v_dot];
    
    A_trans = zeros(6,6);
    A_trans(1:3, 4:6) = R_bi;
    A_trans(4:6, 4:6) = [ 0,  r, -q;
                   -r,  0,  p;
                    q, -p,  0];
end

function [y_hat_trans, H_trans] = TranslationFilterMeasurement(x)
    
    p_d_est = x(3);
    u_est = x(4);
    v_est = x(5);
    
    y_hat_trans = [u_est; v_est; p_d_est];
    
    H_trans = zeros(3, 6);
    H_trans(1, 4) = 1;
    H_trans(2, 5) = 1;
    H_trans(3, 3) = 1;
end
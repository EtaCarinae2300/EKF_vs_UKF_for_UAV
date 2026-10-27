# EKF vs. UKF for UAV State Estimation

This repository contains our **ASEN 5128 Small UAS** project comparing the **Extended Kalman Filter (EKF)** and **Unscented Kalman Filter (UKF)** for nonlinear fixed-wing UAV state estimation in **GPS-denied environments**. The MATLAB simulation includes a 6-DoF aircraft model, autopilot, orbit/straight-line guidance, simulated inertial and non-GPS sensors, state-estimation algorithms, and flight visualization.

> **Current recommended configuration:** Use the **true aircraft state for guidance and control**. The estimated state is still under development and should **not** be used as autopilot feedback. Sensor-tuning changes below are recommendations for a cleaner simulation baseline; they are not automatically applied by this README.

## Navigation

- [Dependencies](#dependencies)
- [How to Execute Code](#how-to-execute-code)
- [True-State Control (Important)](#true-state-control-important)
- [Flight Navigation and Guidance](#flight-navigation-and-guidance)
- [Sensor Tuning to Reduce Outliers](#sensor-tuning-to-reduce-outliers)
- [Folder Structure](#folder-structure)
- [Major Files of Relevance](#major-files-of-relevance)
- [Current Limitations](#current-limitations)

## Dependencies

- **MATLAB**
- **Optimization Toolbox** — `CalculateTrim.m` uses `fmincon` to obtain a trimmed aircraft state.
- Functions used by the filters/tests must also be available in your MATLAB installation, including `wrapToPi` and (for `EKF_TMT.m`) `chi2inv` from Statistics and Machine Learning Toolbox.

No Python environment or external dataset is required for the main MATLAB simulation.

## How to Execute Code

1. Open MATLAB and set the **current folder** to this repository's root directory.
2. Run `startup.m` **before** running other scripts. It adds the project subfolders to the MATLAB path.
3. Run `Run_Project.m` to simulate the T-Twistor aircraft with guidance, autopilot, sensors, and the selected estimator.

```matlab
startup
Run_Project
```

### Main simulation settings

Modify the flags near the top of [`Run_Project.m`](Run_Project.m):

```matlab
ANIMATE_FLAG = 0;         % Set to 1 to animate the aircraft
CONTROL_FLAG = SLC;      % SLC or FEED autopilot
ESTIM_FLAG = UKF;        % UKF selects the hybrid UKF/EKF implementation
ESTIM_CONTROL_FLAG = 0;  % IMPORTANT: 0 = TRUE STATE for control
```

The default simulation duration is `Tfinal = 1000` seconds. Reduce it (for example, to `100` seconds) for a faster initial check. The main script generates trajectory/state plots and prints state-estimation RMSE values; treat estimator metrics as **experimental**, not validated performance claims.

**Other entry points:** [`EKF_TMT.m`](EKF_TMT.m) is an experimental Monte Carlo EKF consistency test using NEES and NIS. It needs further work before use; see [Current Limitations](#current-limitations).

## True-State Control (Important)

In [`Run_Project.m`](Run_Project.m), keep:

```matlab
ESTIM_CONTROL_FLAG = 0;
```

With this setting, the autopilot and guidance logic use the simulated **true state** (`aircraft_array`), rather than `aircraft_state_est`. Do **not** set `ESTIM_CONTROL_FLAG = 1` until the estimators have been debugged and validated.

**Note:** Setting this flag to `0` does **not** disable estimator calculations or estimator comparison plots. It only prevents the estimated state from driving the aircraft. The code currently still executes the selected filter, and unstable estimates may affect its diagnostics or execution.

## Flight Navigation and Guidance

The default flight mission is **orbit following** through [`OrbitGuidance.m`](Guidance%20Algos/OrbitGuidance.m). In [`Run_Project.m`](Run_Project.m), the parameters are:

```matlab
gvf_speed = 18;                  % Desired airspeed [m/s]
gvf_radius = 1000;               % Orbit radius [m]
gvf_center = [5000;5000;-1805]; % Orbit center [North; East; Down] [m]
gvf_flag = 1;                    % Orbit direction
gvf_gains.kr = 0.01;             % Radial guidance gain
gvf_gains.kz = 0.001;            % Vertical guidance gain
```

To follow a different orbit, adjust its center, radius, speed, and direction. The guidance output supplies altitude, climb-rate, course, course-rate feedforward, and airspeed commands to the autopilot.

A **straight-line path-following** option is also available in [`StraightLineGuidance.m`](Guidance%20Algos/StraightLineGuidance.m). To use it, define a point and direction for the desired line, then replace the `OrbitGuidance(...)` call inside the simulation loop with an appropriately configured `StraightLineGuidance(pos_line, dir_line, aircraft_state_con(1:3,i), kpath, chi_inf, gvf_speed)` call. Use a nonzero horizontal component for `dir_line`, as the straight-line guidance function divides by its horizontal magnitude.

For the recommended true-state configuration, `aircraft_state_con` is populated from the true aircraft state.

## Sensor Tuning to Reduce Outliers

The current [`OpticalFlowSensor.m`](Params/OpticalFlowSensor.m) intentionally uses heavy-tailed noise, has a **1% chance of replacing both velocity measurements with zero**, and scales the noise with **absolute altitude**. At the default flight altitude, that can produce extremely large velocity-measurement errors. A cleaner baseline is recommended while debugging estimation.

### Recommended changes

| File | Current behavior / setting | Recommended change | Purpose |
| --- | --- | --- | --- |
| [`SensorParametersTtwistor.m`](Params/SensorParametersTtwistor.m) | `sig_pixel` and `sig_mag` are overridden in `Run_Project.m` | Define `sig_pixel = 0.005`, `sig_mag = 100`, and `optical_flow_dropout_prob = 0` in the central sensor-parameter file | Use a consistent, lower-noise synthetic baseline |
| [`Run_Project.m`](Run_Project.m) | `sensor_params.sig_pixel = 0.05` and `sensor_params.sig_mag = 500` | Remove those overrides after moving settings to `SensorParametersTtwistor.m` | Prevent the main script from undoing the central tuning |
| [`OpticalFlowSensor.m`](Params/OpticalFlowSensor.m) | Heavy-tailed noise and hard-coded 1% dropout | Use Gaussian `randn` noise and a configurable dropout probability (start at `0`) | Eliminate intentional extreme outliers during baseline testing |
| [`OpticalFlowSensor.m`](Params/OpticalFlowSensor.m) | Noise multiplied by `-aircraft_state(3)` (absolute altitude) | Use altitude **above ground**: `max(1, -aircraft_state(3) - sensor_params.h_ground)` | Scale optical-flow uncertainty to the distance from the ground |
| [`EKF_Filtering.m`](Estimation/EKF_Filtering.m) and [`UKF_EKF_Hybrid_Filtering.m`](Estimation/UKF_EKF_Hybrid_Filtering.m) | Optical-flow measurement covariance uses absolute altitude | Use estimated **height above ground** in `R_xy` | Keep the filter's assumed measurement noise consistent with the sensor model |

Add these settings to [`SensorParametersTtwistor.m`](Params/SensorParametersTtwistor.m):

```matlab
%%% Optical flow / magnetometer: clean simulation baseline
sensor_params.sig_pixel = 0.005;
sensor_params.sig_mag = 100;                 % nT
sensor_params.optical_flow_dropout_prob = 0; % 0 = disable synthetic dropouts
```

For the optical-flow generator, the corresponding logic can be:

```matlab
h_agl = max(1, -aircraft_state(3) - sensor_params.h_ground);
sigma_uv = sensor_params.sig_pixel * h_agl;
u_meas = aircraft_state(7) + sigma_uv * randn;
v_meas = aircraft_state(8) + sigma_uv * randn;

if rand < sensor_params.optical_flow_dropout_prob
    u_meas = 0;
    v_meas = 0;
end
```

For the translational EKF / hybrid UKF-EKF measurement covariance, replace the existing `R_xy` expression with:

```matlab
h_agl_est = max(1, h_hat - sensor_params.h_ground);
R_xy = (sensor_params.sig_pixel * h_agl_est)^2;
```

**Why these values?** With a simulated altitude of about **1805 m** and `h_ground = 1655 m`, height above ground is approximately **150 m**. The suggested `sig_pixel = 0.005` gives a nominal optical-flow velocity-noise standard deviation of **0.75 m/s**, rather than the current much larger absolute-altitude-based scale. `sig_mag = 100 nT` is a lower-noise **synthetic test setting**, not a sensor-hardware specification. These are starting points only; validate using repeatable simulation runs before making claims about EKF/UKF accuracy.

**Keep** the accelerometer, gyro, and pressure noise/bias settings unchanged initially to isolate the effect of the optical-flow fix. When moving beyond this clean baseline, reintroduce realistic sensor dropout/outlier behavior and use robust measurement rejection or downweighting in the filters rather than assuming all measurements are Gaussian.

### Applying the supplied tuning patch

A separate `sensor_tuning.patch` accompanies this README. If you want to make the recommended changes in your local checkout, place the patch in the repository root and run:

```bash
git apply --check sensor_tuning.patch
git apply sensor_tuning.patch
```

The patch updates the two sensor-related files, removes the overriding values in `Run_Project.m`, and aligns the optical-flow measurement covariance in the EKF, hybrid UKF/EKF, and standalone UKF. **It has been checked for patch applicability against the uploaded text, but the MATLAB simulation has not been executed or numerically validated.**

## Folder Structure

| Folder | Purpose |
| --- | --- |
| [`Control Law/`](Control%20Law/) | Simple successive-loop-closure autopilot, feedforward autopilot, and control-gain calculations. |
| [`Dynamics/`](Dynamics/) | Nonlinear aircraft equations of motion, aerodynamic forces/moments, inertia calculations, and integration. |
| [`Estimation/`](Estimation/) | EKF, hybrid UKF/EKF, standalone UKF, and reference state estimators. |
| [`Guidance Algos/`](Guidance%20Algos/) | Orbit and straight-line path guidance/navigation algorithms. |
| [`Kinematics/`](Kinematics/) | Coordinate transformations, rotation matrices, wind angles, and flight-path calculations. |
| [`Params/`](Params/) | T-Twistor aircraft data, sensor parameters, simulated GPS/IMU/magnetometer/optical-flow measurements, and saved control gains. |
| [`Plot Visualization/`](Plot%20Visualization/) | Trajectory/state plotting and aircraft animation. |
| [`Trim Calcs/`](Trim%20Calcs/) | Straight-flight and coordinated-turn trim calculations. |

## Major Files of Relevance

- [`startup.m`](startup.m) — Adds the required subfolders to the MATLAB search path.
- [`Run_Project.m`](Run_Project.m) — Main simulation loop, guidance, autopilot, sensors, estimator selection, plotting, and RMSE calculations.
- [`EKF_TMT.m`](EKF_TMT.m) — Experimental Monte Carlo NEES/NIS consistency-test script.
- [`Estimation/EKF_Filtering.m`](Estimation/EKF_Filtering.m) — EKF-based attitude and translational estimation.
- [`Estimation/UKF_EKF_Hybrid_Filtering.m`](Estimation/UKF_EKF_Hybrid_Filtering.m) — Hybrid estimator used when `ESTIM_FLAG = UKF` in the main script.
- [`Estimation/UKF_Filtering.m`](Estimation/UKF_Filtering.m) — Separate full-state UKF implementation; not the filter selected by the default main-script `UKF` branch.
- [`Params/SensorParametersTtwistor.m`](Params/SensorParametersTtwistor.m) — Sampling rates and base sensor-noise parameters.
- [`Params/OpticalFlowSensor.m`](Params/OpticalFlowSensor.m) — Simulated body-frame velocity measurements and outlier behavior.
- [`Params/InertialSensors.m`](Params/InertialSensors.m) — Simulated accelerometer, gyroscope, and pressure measurements.
- [`Params/MagnetometerSensor.m`](Params/MagnetometerSensor.m) — Simulated body-frame magnetic-field measurements.
- [`Dynamics/AircraftEOM.m`](Dynamics/AircraftEOM.m) — Nonlinear 12-state aircraft equations of motion.
- [`Guidance Algos/OrbitGuidance.m`](Guidance%20Algos/OrbitGuidance.m) — Default orbit-navigation commands.
- [`Guidance Algos/StraightLineGuidance.m`](Guidance%20Algos/StraightLineGuidance.m) — Alternate straight-line-navigation commands.
- [`Plot Visualization/PlotSimulationWithCommands.m`](Plot%20Visualization/PlotSimulationWithCommands.m) — Visualizes simulated aircraft states and commanded values.
- [`Trim Calcs/CalculateTrim.m`](Trim%20Calcs/CalculateTrim.m) — Computes straight-flight trim using `fmincon`.

## Current Limitations

- **Estimated-state feedback is not validated.** Keep `ESTIM_CONTROL_FLAG = 0` for control using true aircraft state.
- **The EKF/UKF comparison remains experimental.** The proposed noise adjustments may make the simulated measurements cleaner but do not, by themselves, resolve estimator-model or filter-consistency issues.
- **`EKF_TMT.m` needs repair before interpreting its NEES/NIS results:** It requests a third covariance output from `EKF_Filtering`, while the current filter signature declares only two outputs. The test also contains placeholder/approximation comments.
- The main script has GPS generation commented out; the GPS-based `SIMPLE` and `SMOOTH` estimator branches should not be assumed to work without restoring the appropriate GPS input.

---

*Project: ASEN 5128 — Small Unmanned Aircraft Systems | MATLAB simulation and nonlinear state estimation.*

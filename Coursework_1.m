%% Define and initialisation Constants
Define_Constants;

tau = 0.5;
S_a = 5;
S_cphi_a = 0.01;
S_cf_a = 0.04;
% GNSS pseudo-range std per axis (m)
sigma_pos = 10; 
% GNSS pseudo-range rate std per axis (m/s)
sigma_vel = 0.05; 
% clock offset std (m)
sigma_p = 100000; 
% clock drift std (m/s)
sigma_r = 200; 


% position uncertainty in each direction
sigma_r_DR = 10; 
% velocity uncertainty in each direction
sigma_v_DR = 0.1; 
% DR velocity error power spectral density (PSD)
S_DR = 0.2;



%% Read in data files
Pseudo_ranges = table2array(readtable("Pseudo_ranges.csv"));
Pseudo_range_rates = table2array(readtable("Pseudo_range_rates.csv"));
pseudo_ranges = Pseudo_ranges(2:end, 2:end);
pseudo_range_rates = Pseudo_range_rates(2:end, 2:end);

time = Pseudo_ranges(2:end, 1);
sat_ids = Pseudo_ranges(1, 2:end);

T = length(pseudo_range_rates);


%% Step 1 and 2: Least squares estimation for initial position, velocity, 
% clock offset and drift
[ls_pos_ECEF, ls_vel_ECEF, ls_clock_offs, ls_clock_drift, outlier_idx] = ...
    ls_est(pseudo_ranges, pseudo_range_rates, time, sat_ids);


%% Step 3: Outlier removal and run LS again if necessary
%  Removes outlier satellite if detected in previous step
if outlier_idx > 0
    pseudo_ranges(:,outlier_idx) = [];
    pseudo_range_rates(:,outlier_idx) = [];
    sat_ids(outlier_idx) = [];

    [ls_pos_ECEF, ls_vel_ECEF, ls_clock_offs, ls_clock_drift] = ls_est(...
        pseudo_ranges, pseudo_range_rates, time, sat_ids);
end

N_sat = length(sat_ids);


%% Step 4: Initialise the GNSS Kalman filter estimates
%  This uses the LS solution as the starting estimate for receiver 
%  position, clock offset and clock drift

% Use LS estimate as starting point
x_est = [ls_pos_ECEF(1,:).'; ls_vel_ECEF(1,:).'; ...
    ls_clock_offs(1); ls_clock_drift(1)];

% Initialise error covariance matrix
P_matrix =  diag([sigma_pos, sigma_pos, sigma_pos, sigma_vel, ...
    sigma_vel, sigma_vel, sigma_p, sigma_r]) ^ 2;


%% Step 5: Compute the GNSS Kalman filter solution
% Uses the initial estimates stated in Step 4
[GNSS_pos, GNSS_vel] = Kalman_GNSS(x_est, P_matrix, ...
    pseudo_ranges, pseudo_range_rates, time, sat_ids, sigma_pos, ...
    sigma_vel, tau);


%% Step 6: Gyro-magnetic integration

% Read in the data
Dead_reckoning = table2array(readtable("Dead_reckoning.csv"));
speed = mean(Dead_reckoning(:, 2:5), 2);
gyro_heading = Dead_reckoning(:, 6);
magnetic_heading = Dead_reckoning(:, 7);
time = Dead_reckoning(:, 1);

% Compute the Gyro-magnetic integrastion solution for the heading
heading_integration = Gyro_Magnetometer(time, ...
    magnetic_heading, gyro_heading);


%% Step 7: DR solution w. gyro-magnetic heading
%  Uses the average wheel speed combined with the heading solution from
%  the gyro-magnetometer integration
[ls_lat_NED, ls_lon_NED, ls_elev_NED, ~] = pv_ECEF_to_NED(ls_pos_ECEF(1,:).', ls_vel_ECEF(1,:).');
lat_init = ls_lat_NED * rad_to_deg;
lon_init = ls_lon_NED * rad_to_deg;
height_const = ls_elev_NED; % Assumes a constant height

% Compute the simple DR solution for horizontal position and velocity
[DR_pos, DR_vel] = simple_DR(lat_init, lon_init, speed, heading_integration, time, height_const);


%% Step 8: DR/GNSS Integrated Navigation Kalman filter solution
%  Uses a Kalman filter to estimate the DR error at each time step using
%  GNSS measurements. 

[DR_GNSS_pos, DR_GNSS_vel] = Kalman_DR_GNSS(GNSS_pos, GNSS_vel, DR_pos, ...
    DR_vel, tau, sigma_r_DR, sigma_v_DR, sigma_pos, sigma_vel, S_DR);


%% Motion Profile File 
% Store data to Motion_Profile file
% 1) Time in seconds
% 2) Geodetic latitude in degrees
% 3) Geodetic longitude in degrees
% 4) North velocity in metres per second
% 5) East velocity in metres per second
% 6) Heading in degrees

motion_profile = [time, DR_GNSS_pos(:,1) * rad_to_deg, ...
    DR_GNSS_pos(:,2) * rad_to_deg, DR_GNSS_vel(:,1), DR_GNSS_vel(:,2), ...
    heading_integration.'];
writematrix(motion_profile, 'output/MotionProfile.csv');


%% Plotting
% We plot the traces as .eps files for higher quality presentation in
% LaTeX

plot_traces(time, GNSS_pos, GNSS_vel, DR_pos, ...
    DR_vel, DR_GNSS_pos, DR_GNSS_vel, heading_integration)






function [GNSS_pos_sol, GNSS_vel_sol] = Kalman_GNSS(x_plus, ...
    P_plus, pseudo_ranges, pseudo_range_rates, time, sat_ids, ...
    sigma_pos, sigma_vel, tau)
% Computes the GNSS Kalman filter navigation solution.
%
% Inputs:
%   x_plus               The initial state estimates. First three elements 
%                        are the ECEF position coordinates, next three are 
%                        the ECEF velocities, the 7th is the clock offset 
%                        estimate and the 8th is the clock drift estimate 
%                        (8x1 column vector)
%   P_plus               The initial error covariance matrix (8x8 matrix)
%   pseudo_ranges        Pseudo-range GNSS measurements to receiver in
%                        meters (TxN_sat matrix)
%   pseudo_range_rates   Pseudo-range rate GNSS measurements to receiver in
%                        meters per second (TxN_sat matrix)
%   time                 Time at each epoch in seconds (Tx1 column vector)
%   sat_ids              ID of each satellite (1xN_sat row vector)
%   sigma_pos            Uncertainty of the initial position per axis in 
%                        meters (scalar)
%   sigma_vel            Uncertainty of the initial velocity per axis in 
%                        meters (scalar)
%   tau                  Timelength of each epoch in seconds (scalar)
%
% Outputs:
%   GNSS_pos_sol         NED latitude and longitude solution in degrees as 
%                        well as elevation in meters (Tx3 matrix)
%   GNSS_vel_sol         Velocity solution in m/s (Tx3 matrix)
Define_Constants;

S_a = 5;
S_cphi_a = 0.01;
S_cf_a = 0.04;

N_sat = length(sat_ids);
T = length(pseudo_ranges);

% Initialise return vectors
GNSS_pos_sol = zeros(T, 3);
GNSS_vel_sol = zeros(T, 3);

for k=1:T
    % Formulate the transition matrix
    Phi_k_1 = eye(8);
    Phi_k_1(1:3,4:6) = tau * eye(3);
    Phi_k_1(7,8) = tau;
    
    % Formulate the system noise covariance matrix
    Q_k_1 = eye(8);
    Q_k_1(1:6,1:6) = [1/3 * S_a * tau^3 * eye(3), 1/2 * S_a * tau^2 * eye(3);...
        1/2 * S_a * tau^2 * eye(3), S_a * tau * eye(3)];
    Q_k_1(7:8,7:8) = [S_cphi_a * tau + 1/3 * S_cf_a * tau^3, 1/2* S_cf_a * tau^2;...
        1/2 * S_cf_a * tau^2, S_cf_a * tau];
    
    % Propagate the state estimates
    x_minus = Phi_k_1 * x_plus;

    % Propagate the error covariance matrix
    P_k_minus = Phi_k_1 * P_plus * Phi_k_1.' + Q_k_1;
    
    % Compute the range from the estimated receiver position to each
    % satellite
    r_hat_a_minus = zeros(N_sat, 1);
    for it=1:2
        for j=1:N_sat
            C_e_I = eye(3) - Omega_ie * r_hat_a_minus(j) / c;
            [sat_pos, ~] = Satellite_position_and_velocity(time(k), sat_ids(j));
            r_hat_a_minus(j) = norm(C_e_I * sat_pos.' - x_minus(1:3));
        end
    end
    
    % Compute the line-of-sight unit vector from the estimated receiver
    % position to each satellite
    u_a_e = zeros(3, N_sat);
    for j=1:N_sat
        C_e_I = eye(3) - Omega_ie * r_hat_a_minus(j) / c;
        [sat_pos, ~] = Satellite_position_and_velocity(time(k), sat_ids(j));
        u_a_e(:, j) = (C_e_I * sat_pos.' - x_minus(1:3)) / r_hat_a_minus(j);
    end
    
    % Compute the range rate from the estimated receiver position to each
    % satellite
    r_dot_hat_a_minus = zeros(N_sat, 1);
    for j=1:N_sat
        C_e_I = eye(3) - Omega_ie * r_hat_a_minus(j) / c;
        [sat_pos, sat_vel] = Satellite_position_and_velocity(time(k), sat_ids(j));
        r_dot_hat_a_minus(j) = u_a_e(:, j).' * (C_e_I * (sat_vel.' + Omega_ie * sat_pos.') - (x_minus(4:6) + Omega_ie * x_minus(1:3)));
    end
    
    % Formulate the measurement matrix
    H_k = [-u_a_e.', zeros(width(u_a_e), 3), ones(width(u_a_e), 1), zeros(width(u_a_e), 1);...
        zeros(width(u_a_e), 3), -u_a_e.', zeros(width(u_a_e), 1), ones(width(u_a_e), 1)];
    
    % Formulate the measurement noise covariance matrix
    R_k = diag([sigma_pos^2 + zeros(1, N_sat), sigma_vel^2 + zeros(1, N_sat)]);
    
    % Compute the Kalman gain matrix
    K_k = P_k_minus * H_k.' / (H_k * P_k_minus * H_k.' + R_k);
    
    % Formulate the measurement innovation vector
    rho = pseudo_ranges(k, :).';
    rho_dot = pseudo_range_rates(k, :).';
    clock_offs_est = x_minus(7);
    clock_drift_est = x_minus(8);
    dz_min = [rho; rho_dot] ...
        - [r_hat_a_minus; r_dot_hat_a_minus] ...
        - [repelem(clock_offs_est, N_sat).'; repelem(clock_drift_est, N_sat).'];
    
    % Compute the new state estimates
    x_plus = x_minus + K_k * dz_min;
    % Compute the new error covariance matrix
    P_plus = (eye(height(K_k)) - K_k * H_k) * P_k_minus;
    
    % Map to degrees and update the return vectors
    [lat, lon, elev, vel] = pv_ECEF_to_NED(x_plus(1:3), x_plus(4:6));
    lat = lat * rad_to_deg;
    lon = lon * rad_to_deg;

    GNSS_pos_sol(k, :) = [lat, lon, elev];
    GNSS_vel_sol(k, :) = vel;
end

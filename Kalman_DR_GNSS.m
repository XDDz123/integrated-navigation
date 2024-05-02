function [DR_GNSS_pos, DR_GNSS_vel] = Kalman_DR_GNSS(GNSS_pos, ...
    GNSS_vel, DR_pos, DR_vel, tau_s, sigma_pos, sigma_vel, ...
    sigma_Gr, sigma_Gv, S_DR)
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the DR/GNSS Integrated Navigation solution using a Kalman filter
% by estimating the Dead Reckoning error at each epoch with the GNSS
% solution.
%
% Inputs:
%   GNSS_pos      GNSS Kalman filter NED position solution in degrees
%                 (Tx3 matrix)
%   GNSS_vel      GNSS Kalman filter velocity solution in m/s (Tx2 matrix)
%   DR_pos        Dead reckoning NED position solution in radians 
%                 (Tx2 matrix)
%   DR_vel        Dead reckoning NED velocity solution in m/s (Tx2 matrix)
%   tau_s         Timelength of each epoch in seconds (scalar)
%   sigma_pos     Uncertainty of the initial position per axis in meters
%                 (scalar)
%   sigma_vel     Uncertainty of the initial velocity per axis in meters
%                 (scalar)
%   sigma_Gr      Uncertainty of the GNSS position measurements per axis 
%                 in meters (scalar)
%   sigma_Gv      Uncertainty of the GNSS velocity measurements per axis 
%                 in m/s (scalar)
%   S_DR          Dead reckoning velocity error power spectral density
%                 (PSD) in m^2s^-3 (scalar)
%
% Outputs:
%   DR_GNSS_pos   Latitude and longitude solution in radians (Tx2 matrix)
%   DR_GNSS_vel   Horizontal velocity solution in m/s (Tx2 matrix)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Define_Constants;

T = length(GNSS_pos);

% Dead reckoning solution (lat/lon in radians)
L_D = DR_pos(:, 1);
lambda_D = DR_pos(:, 2);
v_N_D = DR_vel(:, 1);
v_E_D = DR_vel(:, 2);

% GNSS solution (lat/lon in radians)
L_G = GNSS_pos(:, 1) * deg_to_rad;
lambda_G = GNSS_pos(:, 2) * deg_to_rad;
v_N_G = GNSS_vel(:, 1);
v_E_G = GNSS_vel(:, 2);

L_C = zeros(T, 1); % Corrected latitude solution (rad)
L_C(1) = L_G(1);
lambda_C = zeros(T, 1); % Corrected longitude solution (rad)
lambda_C(1) = lambda_G(1);
v_N_C = zeros(T, 1); % Corrected north vel solution
v_N_C(1) = GNSS_vel(1, 1);
v_E_C = zeros(T, 1); % Corrected east vel solution
v_E_C(1) = GNSS_vel(1, 2);

[R_N, R_E] = Radii_of_curvature(L_C(1));
h = GNSS_pos(:, 3);

x_plus = [0;0;0;0]; % dv_N (north vel error), dv_E (east vel error), dL (lat error), dlambda (lon error)

P_plus = diag([sigma_vel^2, sigma_vel^2, sigma_pos^2 / (R_N + h(1))^2, ...
    sigma_pos^2 / ((R_E + h(1)) * cos(L_C(1))^2)]);

for k=2:T
    [R_N, R_E] = Radii_of_curvature(L_C(k-1));
    
    % Transition matrix
    Phi = eye(4);
    Phi(3, 1) = tau_s / (R_N + h(k-1));
    Phi(4, 2) = tau_s / ((R_E + h(k-1)) * cos(L_C(k-1)));
    
    % System noise covariance matrix
    Q = [S_DR * tau_s, 0, 1/2 * S_DR * tau_s^2 / (R_N + h(k-1)), 0;...
        0, S_DR * tau_s, 0, 1/2 * S_DR * tau_s^2 / ((R_E + h(k-1)) * cos(L_G(k-1)));...
        1/2 * S_DR * tau_s^2 / (R_N + h(k-1)), 0, 1/3 * S_DR * tau_s^3 / (R_N + h(k-1))^2, 0;...
        0, 1/2 * S_DR * tau_s^2 / ((R_E + h(k-1)) * cos(L_G(k-1))), 0, 1/3 * S_DR * tau_s^3 / ((R_E + h(k-1))^2 * cos(L_G(k-1))^2)];
    
    % Propagate state estimates
    x_minus = Phi * x_plus;
    
    % Propagate error covariance matrix
    P_minus = Phi * P_plus * Phi.' + Q;
    
    % Measurement matrix 
    % (assumes pos measurements precede vel measurements)
    H_k = diag([-1, -1], 2) + diag([-1, -1], -2);
    
    % Measurement noise covariance matrix
    R_k = diag([sigma_Gr^2 / (R_N + h(k))^2, sigma_Gr^2 / ((R_E + h(k))^2 * cos(L_G(k))^2), sigma_Gv^2, sigma_Gv^2]);
    
    % Kalman gain matrix
    K_k = P_minus * H_k.' / (H_k * P_minus * H_k.' + R_k);
    
    % Measurement innovation vector
    dz_k_minus = [L_G(k); lambda_G(k); v_N_G(k); v_E_G(k)] ...
        - [L_D(k); lambda_D(k); v_N_D(k); v_E_D(k)] ...
        - H_k * x_minus;
    
    % Update state estimates
    x_plus = x_minus + K_k * dz_k_minus;
    
    % Update error covariance matrix
    P_plus = (eye(length(K_k)) - K_k * H_k) * P_minus;
    
    % Add to corrected solution
    v_N_C(k) = v_N_D(k) - x_plus(1);
    v_E_C(k) = v_E_D(k) - x_plus(2);
    L_C(k) = L_D(k) - x_plus(3);
    lambda_C(k) = lambda_D(k) - x_plus(4);
end

DR_GNSS_pos = [L_C, lambda_C];
DR_GNSS_vel = [v_N_C, v_E_C];

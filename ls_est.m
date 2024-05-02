function [ls_pos_ECEF, ls_vel_ECEF, ls_clock_offs, ls_clock_drift, outlier_idx] = ls_est(Pseudo_ranges, Pseudo_range_rates, time, sat_ids)
% Computes the least squares solution for the position, velocity, 
% clock offset and drift at each epoch using GNSS measurements
%
% Inputs:
%   Pseudo_ranges        Pseudo-range GNSS measurements to receiver in
%                        meters (TxN_sat matrix)
%   Pseudo_range_rates   Pseudo-range rate GNSS measurements to receiver in
%                        meters per second (TxN_sat matrix)
%   time                 Time at each epoch in seconds (Tx1 column vector)          
%   sat_ids              ID of each satellite (1xN_sat row vector)
%
% Outputs:
%   ls_pos_ECEF          ECEF least squares position solution (Tx3 matrix)
%   ls_vel_ECEF          ECEF least squares velocity solution (Tx3 matrix)
%   ls_clock_offs        Least squares clock offset solution (Tx1 column
%                        vector)
%   ls_clock_drift       Least squares clock drift solution (Tx1 column
%                        vector)
%   outlier_idx          The index of an outlier if one is detected
%                        (otherwise -1)
Define_Constants;

N_sat = length(sat_ids);
ls_pos_ECEF = zeros(length(Pseudo_ranges), 3);
ls_vel_ECEF = zeros(length(Pseudo_range_rates), 3);
ls_clock_offs = zeros(length(Pseudo_range_rates), 1);
ls_clock_drift = zeros(length(Pseudo_range_rates), 1);

x_plus = [0; 0; 0; 0];
cell_pos = [0; 0; 0]; % Initial prediction
x_plus_rate = [0; 0; 0; 0];
v_user_pred = [0; 0; 0]; % Initial prediction
outlier_idx = -1; % Becomes positive if an outlier is detected

for k = 1 : length(Pseudo_ranges)

    % Iterate until convergence
    while 1
        last_x_plus = x_plus;
        sat_pos = zeros(N_sat, 3);
        sat_vel = zeros(N_sat, 3);
        
        % Compute the position and velocity of the satellites 
        for i=1:N_sat
            [sat_pos(i,:), sat_vel(i,:)] = ...
                Satellite_position_and_velocity(time(k), sat_ids(i));
        end
        
        % Compute the ranges from the approximate user position to each 
        % satellite, taking into account the Sagnac effect
        iterations = 2;
        ranges = zeros(1, N_sat);
        range_rates = zeros(1, N_sat);
        
        for it=1:iterations
            for i=1:N_sat
                C_e_I = eye(3) - Omega_ie * ranges(i) / c;
                vec = C_e_I * transpose(sat_pos(i, :)) - cell_pos;
                ranges(i) = sqrt(dot(vec, vec));
            end
        end
        
        % Compute the line-of-sight unit vectors to each satellite
        line_of_sights = zeros(N_sat, 3);
        for i=1:N_sat
            C_e_I = eye(3) - Omega_ie * ranges(i) / c;
            line_of_sights(i, :) = (C_e_I * transpose(sat_pos(i, :)) ...
                - cell_pos) / ranges(i);
        end
        
        % Compute the range rates from the user position to each satellite
        for i=1:N_sat
            C_e_I = eye(3) - Omega_ie * ranges(i) / c;
            range_rates(i) = line_of_sights(i, :) * (C_e_I * ...
                (transpose(sat_vel(i, :)) + Omega_ie * ...
                transpose(sat_pos(i, :))) - ...
                (v_user_pred + Omega_ie * cell_pos));
        end
        
        % Formulate the predicted range and clock offset state vector
        clock_offs_pred = x_plus(4);
        x_pred = [cell_pos; clock_offs_pred];
        
        % Formulate the range measurement innovation vector
        dz_min = zeros(N_sat, 1);
        for i=1:N_sat
            dz_min(i) = Pseudo_ranges(k, i) - ranges(i) - clock_offs_pred;
        end
        
        % Formulate the range measurement matrix
        H_G_e = zeros(N_sat, 4);
        for i=1:N_sat
            H_G_e(i, :) = [-line_of_sights(i, :) 1];
        end
        
        % Formulate the predicted range rate and clock drift state vector
        clock_drift = x_plus_rate(4);
        x_pred_rate = [v_user_pred; clock_drift];

        % Formulate the rate measurement innovation vector
        dz_min_rate = zeros(N_sat, 1);
        for i=1:N_sat
            dz_min_rate(i) = Pseudo_range_rates(k, i) - range_rates(i) ...
                - clock_drift;
        end

        % Compute the new estimates using unweighted least squares
        x_plus = x_pred + pinv(H_G_e) * dz_min;
        x_plus_rate = x_pred_rate + pinv(H_G_e) * dz_min_rate;
        
        % Populate the return vectors
        ls_pos_ECEF(k, :) = x_plus(1:3);
        ls_vel_ECEF(k, :) = x_plus_rate(1:3);
        ls_clock_offs(k) = x_plus(4);
        ls_clock_drift(k) = x_plus_rate(4);
        cell_pos = x_plus(1:3);
        v_user_pred = x_plus_rate(1:3);
        

        % Convergence check
        if norm(x_plus - last_x_plus) < 0.001
            break;
        end

        % Check for outliers
        m = length(sat_ids);
        v = (H_G_e * pinv(H_G_e) - eye(m)) * dz_min;
    
        C_v = (eye(m) - H_G_e * pinv(H_G_e)) * 25;
        
        outlier_detected = any(abs(v) > sqrt(diag(C_v)) * 6);
        if outlier_detected
            [~, outlier_idx] = max(abs(v));
        end
    end
end



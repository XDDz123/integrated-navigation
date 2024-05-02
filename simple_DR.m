function [latlon_DR, v_inst_DR] = simple_DR(lat_init, lon_init, speed, heading, time, height_const)
% Inputs:
%   lat_init      Initial latitude in degrees (scalar)
%   lon_init      Initial longitude in degrees (scalar)
%   speed         Speed measurements in m/s (Tx1 column vector)
%   heading       Heading measurements in degrees (Tx1 column vector)
%   time          Time at each epoch (Tx1 column vector)
%   height_const  Assumed NED height in m (scalar)
%
% Outputs:
%   latlon_DR     Latitude and longitude dead reckoning solution (Tx2
%                 matrix)
%   v_inst_DR     Damped instantaneous horizontal velocity dead reckoning 
%                 solution (Tx2 matrix)
Define_Constants;

h = height_const;
T = length(time);

v_avg_DR = zeros(T, 2);
v_inst_DR = zeros(T, 2);

% Map heading measurements to radians
psi = heading * deg_to_rad;

% Initialise the damped instantaneous velocity matrix
v_inst_DR(1, :) = [
    speed(1) * cos(deg_to_rad * heading(1)),...
    speed(1) * sin(deg_to_rad * heading(1))];

% Initialise the dead reckoning position solution matrix
latlon_DR = zeros(T, 2);
latlon_DR(1,:) = [lat_init, lon_init] * deg_to_rad;

for k=2:T
    % Compute the average velocity in each direction from the speed and
    % heading measurements
    v_avg_DR(k, :) = 1/2 * [cos(psi(k)) + cos(psi(k-1)); ...
        sin(psi(k)) + sin(psi(k-1))] * speed(k); 

    % Compute the dead reckoning position using the estimated average
    % velocity during the epoch
    [R_N, R_E] = Radii_of_curvature(latlon_DR(k-1, 1));
    latlon_DR(k, 1) = latlon_DR(k-1, 1) + v_avg_DR(k, 1) * (time(k) - time(k-1)) / (R_N + h);
    latlon_DR(k, 2) = latlon_DR(k-1, 2) + v_avg_DR(k, 2) * (time(k) - time(k-1)) / ((R_E + h) * cos(latlon_DR(k, 1)));
    
    % Compute the damped instantaneous velocity
    v_inst_DR(k, :) = [
        1.7 * v_avg_DR(k, 1) - 0.7 * v_inst_DR(k-1, 1);...
        1.7 * v_avg_DR(k, 2) - 0.7 * v_inst_DR(k-1, 2)];
end

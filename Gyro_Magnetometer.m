function headingIntegrate = Gyro_Magnetometer(time,heading,gyro_heading)

% Define and initialisation Constants
Define_Constants 
headingResults = zeros(1,size(time,1));
%propagation interval
tau = 0.5; 
%Gyroscope bias standard deviation
gyroscopeBias = 1; 
%Magnetic Compass noise error standard deviation
compassNoiseError = 4*deg_to_rad; 
%Gyroscope random noise PSD
randomNoise = e-4; 
%Gyroscope bias random walk PSD
baisRandomWalk = 3e-6; 
% Kalman Filter state vector
kalmanSilterStateVector = [0;0];
% error covariance matrix
convarianceMatrixError = [compassNoiseError^2 0;0 gyroscopeBias^2];

for epoch=1:size(time,1)

    % Qb Cal. transition matrix
    phik = [1 tau;0 1];

    % Qc Cal. system noise covariance matrix
    Qk1 = [randomNoise*tau+1/3*baisRandomWalk*tau^3 ...
                  1/2*baisRandomWalk*tau^2;
                  1/2*baisRandomWalk*tau^2 baisRandomWalk*tau];
 
    % Qd state estimates
    xk = phik*kalmanSilterStateVector;
    % Qe error covariance matrix
    Pk = phik * convarianceMatrixError * phik' ...
        + Qk1;
    
    % Qf Formulate measurement matrix
    Hk = [-1 0];
    % Qg Cal noise covariance matrix
    Rk = diag(compassNoiseError^2);
    % Qh Cal. Kalman gain matrix
    Kk = Pk*Hk'/(Hk*Pk*Hk' + Rk);
    % Qi Cal. innovation vector
    dz = heading(epoch) - gyro_heading(epoch) - Hk*xk;
    % Qj state estimates update
    xkNew = xk + Kk*dz;
    % Qk error covariance matrix update
    PkNew = (eye(size(Pk,1)) - Kk * Hk) * Pk;
    
    
    % variables update
    kalmanSilterStateVector = xkNew;
    convarianceMatrixError = PkNew;
    headingResults(:,epoch) = gyro_heading(epoch) - xkNew(1);
end

headingIntegrate = headingResults;

end


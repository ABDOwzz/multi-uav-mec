% This function computes the channel gain of the link between the user and the UAV.
% rUser is a 1*3 vector denoting 3-D coordinates of the user.
% rUAV is a 1*3 vector denoting 3-D coordinates of the UAV.

function [g, L] = chanGain(rUAV, rUser)
    fc = 2e9;                                   % Carrier frequency (Hz)
    c = 3e8;                                    % Speed of light (m/sec)
    alpha = 2;                                 % Path loss exponent
    % Consider the urban scenario
    % The coefficients below are based on 
    a = 11.95;
    b = 0.14;
    etaLoS = 10^(1/10);
    etaNLoS = 10^(20/10);
    
    % LoS distance between UAV and user
    d = sqrt((rUAV(1)-rUser(1))^2+(rUAV(2)-rUser(2))^2+(rUAV(3)-rUser(3))^2);
    % Elevation angle of UAV respective to the user.
    theta = 180/pi*asin(abs(rUAV(3)-rUser(3))/d);
    
    % Compute the probability of LoS and NLoS
    pLoS = 1/(1+a*exp(-b*(theta-a)));
    pNLoS = 1-pLoS;
    
    % Average path loss
    L = pLoS*etaLoS*(4*pi*fc*d/c)^alpha+pNLoS*etaNLoS*(4*pi*fc*d/c)^alpha;
    g = 1/L;                                    % Channel gain of the link
end


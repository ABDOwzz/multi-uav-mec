% This function computes the channel gain of the link between the user and the UAV.
% rUser is a 1*3 vector denoting 3-D coordinates of the user.
% rUAV is a 1*3 vector denoting 3-D coordinates of the UAV.

function [g, L] = chanGainFSPL(rUAV, rUser)
    fc = 2e9;                                   % Carrier frequency (Hz)
    c = 3e8;                                    % Speed of light (m/sec)
    alpha = 2;                                 % Path loss exponent
    
    % LoS distance between UAV and user
    d = sqrt((rUAV(1)-rUser(1))^2+(rUAV(2)-rUser(2))^2+(rUAV(3)-rUser(3))^2);
    L = (4*pi*fc*d/c)^alpha;
    g = 1/L;                                    % Channel gain of the link
end


clear all
close all
% This script describes a NOMA-assisted MEC system on UAV according to doc 'Model 1.pages'.
% Created and copyrighted by ZHANG Xiaochen at 10:22 a.m., Apr. 12, 2019.

%% Setup of the Model

% Consider U terrestrial users located in a area
U = 8;                                                 % Number of users/SMDs

% The size of the area is defined by:
lengthArea = 1000;                            % Length of the coverage area (m)
widthArea = 1000;                             % Width of the coverage area (m)

% U SMDs are uniformly distributed in the rectangular coverage area.
rUser = [lengthArea*rand(U, 1),widthArea*rand(U, 1), zeros(U, 1)];
% rUser = [100*(2:2:8).', 100*8*ones(U, 1), zeros(U, 1)];
% rUser = [[400,400,600,600].', [600,800,600,800].',zeros(U, 1)];
% rUser = [[100,500,800,600].', [600,600,600,800].',zeros(U, 1)];
rUser = [kron(100*(2:2:8).',ones(2,1)), ...
    repmat([600, 800].', 4, 1),...
    zeros(U, 1)];



% Divide the time interval into N tme slots.
T = 100;                                               % Length of interval
N = 20;                                                 % Number of time slots
tau = T/N;                                            % Duration of each time slot

% Set the initial trajectory of UAV
M = 2;
hUAV = 50;                                         % Flying altitude of UAV
rI = [0, 0, hUAV; 0, 0, hUAV];              % Intial position of the UAV (M*3)
rF = [lengthArea, 0, hUAV; 0, widthArea, hUAV];               % Final position of the UAV

rUAV = zeros(N, 3, M);                             % Trajectory vector of the UAV

% rUAV(time, dimension, UAV index) is a N*3*M matrix of which each describes the 3-D coordinates of UAVs
% at each time slot.

for m = 1:M
    rUAV(N, :, m) = rF(m, :);
end
for n = 1:N-1
    for m = 1:M
        rUAV(n, :, m) = rI(m, :)+n/N*(rF(m, :)-rI(m, :));
    end
end

% Other parameters involved in the model.

omegaUser = 0.8;                        % Weight factor of users
omegaUAV = 1-omegaUser;         % Weight factor of UAVs

N0 = 1e-17;                                 % Power spectral density of AWGN (W/Hz)
B = 4e6;                                        % System bandwidth (Hz)

D = 120*ones(U, 1);                     % Number of bits for each SMDs to finish (Mbit)
kappa = 1e-28;                            % Effective swiched capacitance of CPU

fmaxUAV = 20;                            % Maximum CPU frequency of MEC server (GHz)
fminUAV = 0;                               % Minimum CPU frequency of SMDs (GHz)
CUAV = 1e3;                                % Processing density of MEC server (Hz/bit)

fmaxUser = 1;                             % Maximum CPU frequency of SMDs (GHz)
fminUser = 0;                              % Minimum CPU frequency of SMDs (GHz)
CUser = 1e3;                               % Processing density of SMD (Hz/bit)

J = log(2)/(B*tau);
J = J*1e6;

vmax = 15;                                   % Maximum flying speed of UAV (m/sec)
alpha = 2;                                     % Path loss factor
fc = 2e9;                                       % Carrier frequency (Hz)
c = 3e8;                                        % Speed of light (m/sec)
Q = (4*pi*fc/c)^alpha;
gamma = 1e-4;
W = 10;                                         % Mass of UAV (kg)

c1 = omegaUser*kappa*tau*1e27;
c2 = omegaUser*B*N0*1e9;
c3 = omegaUAV*kappa*tau*1e27;

%% Approximation of channel gain by curve fitting
a = 11.95;
b = 0.14;

etaLoS = 10^(1/10);
etaNLoS = 10^(20/10);

d = hUAV:1000;
L = (...
    1./(...
    1+a*exp(...
    -180*b/pi*asin(hUAV./d)+a*b...
    )...
    )*(etaLoS-etaNLoS)+etaNLoS...
    )...
    .*(4*pi*fc*d/c).^2;

[curvePoly, goodness, output] = fit(d.', L.', 'poly2');

p1 = curvePoly.p1;
p2 = curvePoly.p2;
p3 = curvePoly.p3;

%% Joint Optimization of Resource Allocation and Trajectory Scheduling

maxIteration = 15;
utilitySeq = zeros(maxIteration, 1);
numIteration = 1;                            % Number of iteration 

pmax = 10;                                      % Maximum offloading power for SMDs.
pmin = 0;
xmin = 0;

%% Find the channel gain for each link between UAVs and users in each time slot
g = zeros(U, N, M);                             % Channel gain matrix
L = zeros(U, N, M);                             % Path loss matrix
% g is a U*N matrix denoting the channel gains for all links throughout the interval.

for i = 1:U
    for n = 1:N
        for m = 1:M
            [g(i, n, m), L(i, n, m)]  = chanGain(rUAV(n, :, m), rUser(i, :));
        end
    end
end

%% Resource Allocation (RA) for Offloading and Computation

cvx_begin
    variable fUser(U, N) nonnegative
    variable fUAV(U, N, M) nonnegative
    variable threshold
    variable userRate(U, N, M)
    expressions offloadPower(U, N M)
    for m = 1:M
        for n = 1:N
            for i = 1:U
                offloadPower(i, n, m) = B*N0/(g(i, n, m)*U)*(exp(U*J*userRate(i, n, m))-1);
            end
        end
    end

    eta = omegaUser*sum(sum(sum(offloadPower, 3)))+c1*sum(sum(pow_p(fUser, 3))) +...
            c3*sum(sum(sum(pow_p(fUAV, 3))));
    minimize( threshold )
    subject to
        % C0: Objective function
        eta<= threshold;
        % C1: CPU frequency constraint on SMDs
        fminUser <= fUser;
        fUser <= fmaxUser;

        % C2: CPU frequency constraint on MEC server
        for m = 1:M
            fminUAV <= sum(fUAV(:, :, m), 1);
            sum(fUAV(:, :, m), 1) <= fmaxUAV;
        end

        % C3: Offloading rate/power constraint on SMDs
        sum(offloadPower, 3) <= pmax;
        userRate>=0;

        % C4: Accomplishment of computation tasks
        D <= 1000*tau/CUser*sum(fUser, 2)+1000*tau/CUAV*sum(sum(fUAV(:, 2:end, :), 3), 2);

        % C5&C6: Causal constraint on offloading and edge computing

        for i = 1:U    
            for n = 2:N
                for m = 1:M
                    sum(userRate(i, (1:n-1), m)) >= tau/CUAV*1000*sum(fUAV(i, (2:n), m));
                end
            end
        end

cvx_end

%% Iteration begins
while 1    
    %% Optimization of UAV Trajectory Scheduling
    chi = zeros(U, N, M);
    % chi is an U*N*M matrix.
    for m= 1:M
        for i = 1:U
            for n = 1:N
                chi(i, n, m) = omegaUser*1e10*B*N0/U*(exp(U*J*userRate(i, n, m))-1);
            end
        end
    end

    for k = 1:U
        for n = 1:N
            for m = 1:M
                if chi(k, n, m) <= 1e-6
                    chi(k, n, m) = 0;
                end
            end
        end
    end

    cvx_begin
        variable r1(N, M)
        variable r2(N, M)
        variable threshold2
        Emove = 0;
        for m = 1:M
            Emove = Emove + square_pos((r1(1, m)-rI(m, 1))/tau)+square_pos((r2(1, m)-rI(m, 2))/tau);
            for n = 2:N
                Emove = Emove+square_pos((r1(n, m)-r1(n-1 ,m))/tau)+square_pos((r2(n, m)-r2(n-1 ,m))/tau);
            end
        end
        Emove = omegaUAV*gamma*tau*W/2*Emove;

        Eoffload = 0;
        for n = 1:N
            for i = 1:U
                for m = 1:M
                Eoffload = Eoffload+chi(k, n, m).*...
                    (...
                    p1*square_pos...
                    (...
                    norm([r1(n, m),r2(n, m),hUAV]-[rUser(i, 1), rUser(i, 2), 0])...
                    )...
                    +p2*norm([r1(n, m),r2(n, m),hUAV]-[rUser(i, 1), rUser(i, 2), 0])...
                    +p3...
                    )/(1e10);
                end
            end
        end


        minimize threshold2
        subject to

            Emove + Eoffload <= threshold2;

            0 <= r1 <= lengthArea;
            0 <= r2 <= widthArea;


            for m = 1:M
            r1(N, m) == rF(m, 1);
            r2(N, m) == rF(m, 2);
                (r1(1, m)-rI(m, 1))^2+(r2(1, m)-rI(m, 2))^2 <= (vmax*tau)^2;
                for n = 2:N
                    (r1(n, m)-r1(n-1, m))^2+(r2(n, m)-r2(n-1, m))^2 <= (vmax*tau)^2;
                end
            end


    cvx_end
    
    % Update the UAV trajectory.
    for m = 1:M
        rUAV(:, 1, m) = r1(:, m);
        rUAV(:, 2, m) = r2(:, m);
    end
    
    %% Find the channel gain for each link between UAVs and users in each time slot
    g = zeros(U, N, M);                             % Channel gain matrix
    L = zeros(U, N, M);                             % Path loss matrix
    % g is a U*N matrix denoting the channel gains for all links throughout the interval.

    for i = 1:U
        for n = 1:N
            for m = 1:M
                [g(i, n, m), L(i, n, m)]  = chanGain(rUAV(n, :, m), rUser(i, :));
            end
        end
    end

    
    %% Resource Allocation (RA) for Offloading and Computation
    cvx_begin
    variable fUser(U, N) nonnegative
    variable fUAV(U, N, M) nonnegative
    variable threshold
    variable userRate(U, N, M)
    expressions offloadPower(U, N M)
    for m = 1:M
        for n = 1:N
            for i = 1:U
                offloadPower(i, n, m) = B*N0/(g(i, n, m)*U)*(exp(U*J*userRate(i, n, m))-1);
            end
        end
    end

    eta = omegaUser*sum(sum(sum(offloadPower, 3)))+c1*sum(sum(pow_p(fUser, 3))) +...
            c3*sum(sum(sum(pow_p(fUAV, 3))));
    minimize( threshold )
    subject to
        % C0: Objective function
        eta<= threshold;
        % C1: CPU frequency constraint on SMDs
        fminUser <= fUser;
        fUser <= fmaxUser;

        % C2: CPU frequency constraint on MEC server
        for m = 1:M
            fminUAV <= sum(fUAV(:, :, m), 1);
            sum(fUAV(:, :, m), 1) <= fmaxUAV;
        end

        % C3: Offloading rate/power constraint on SMDs
        sum(offloadPower, 3) <= pmax;
        userRate>=0;

        % C4: Accomplishment of computation tasks
        D <= 1000*tau/CUser*sum(fUser, 2)+1000*tau/CUAV*sum(sum(fUAV(:, 2:end, :), 3), 2);

        % C5&C6: Causal constraint on offloading and edge computing

        for i = 1:U    
            for n = 2:N
                for m = 1:M
                    sum(userRate(i, (1:n-1), m)) >= tau/CUAV*1000*sum(fUAV(i, (2:n), m));
                end
            end
        end

    cvx_end
    
    %% Update the weighted system 
    % Update the each term of weighted system utility.
    E_user_offload = omegaUser*sum(sum(sum(offloadPower)));
    
    E_user_comp = c1*sum(sum(fUser.^3));
    E_uav_comp = c3*sum(sum(sum(fUAV.^3)));
    E_uav_move = Emove;

    % Update the weighted system utility.
    utilitySeq(numIteration) = E_user_comp+E_user_offload+E_uav_comp+E_uav_move;
    
    % Check  the convergence of the algorithm to break the 'while' loop.
    if numIteration > 1
        if abs(utilitySeq(numIteration)-utilitySeq(numIteration-1)) <= 1
            break
        end
    end
    
    if numIteration == maxIteration
        break
    end
    
    numIteration = numIteration+1
end

% Plot the configuration of the model
figure(1)
scatter(rUser(:, 1), rUser(:, 2));
hold on
stem(rUAV(:, 1, 1), rUAV(:, 2, 1));
hold on
stem(rUAV(:, 1, 2), rUAV(:, 2, 2));

title('Location of SMDs and Trajectory of the UAV')
xlabel('Length/m')
ylabel('Width/m')
axis([0 lengthArea 0 widthArea])

figure(2)
recordOffloadPower = sum(offloadPower, 3);

for i = 1:U
	stem3(i*ones(N, 1), (1:N), recordOffloadPower(i, :));
    hold on
end
xlabel('User Index')
ylabel('Time Slot')
zlabel('Offloading Power')
title('Record of ffloading power of each user')

figure(3)
plot(utilitySeq(1:numIteration))
xlabel('Number Of Iteration')
ylabel('Weighted System Utility')


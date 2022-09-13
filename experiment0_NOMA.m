function [power_NOMA, rUAV_NOMA, utility_NOMA, fUAV_NOMA, fUser_NOMA] = experiment0_NOMA(rUser)

%% Setup of the Model

% Consider U terrestrial users located in a area
U = 8;                                                 % Number of users/SMDs

% The size of the area is defined by:
lengthArea = 1000;                            % Length of the coverage area (m)
widthArea = 1000;                             % Width of the coverage area (m)

% U SMDs are uniformly distributed in the rectangular coverage area.

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
c1 = omegaUser*kappa*tau*1e27;
c2 = omegaUser*B*N0*1e9;
c3 = omegaUAV*kappa*tau*1e27;

vmax = 15;                                   % Maximum flying speed of UAV (m/sec)
alpha = 2;                                     % Path loss factor
fc = 2e9;                                       % Carrier frequency (Hz)
c = 3e8;                                        % Speed of light (m/sec)
Q = (4*pi*fc/c)^alpha;
gamma = 1e-4;
W = 10;                                         % Mass of UAV (kg)

c4 = omegaUser*B*N0;

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
    
    gRank = zeros(U, N, M);
    rank = zeros(U, N, M);
    Psi = zeros(U, N, M);
    beta = zeros(U+2, N, M);
    mu = zeros(U+1, N, M);

    for m = 1:M
        [gRank(:, :, m), rank(:, :, m), Psi(:, :, m)] = chanRank(g(:, :, m), U, N);

        beta(2:(U+1), :, m) = 1./gRank(:, :, m);

        for k = 1:U+1
            mu(k, :, m) = -beta(k+1, :, m)+beta(k, :, m);
        end
    end
    % rank gives a list of user index. Its value denotes the user index with that ranking, 
    % i.e. rank(k, n) means the index of user with the k-th smallest channel gain at time n.
    % Psi helps track the rankings of users, 
    % i.e. Psi(i, n) denotes the ranking of user i at time n.
    % beta is a (U+2)*N matrix.
    mu = mu/(1e9);
    
    % Set initial values for iteration.
    x0 = zeros(U, N, M);
    eta0 = 0;
%% Resource Allocation (RA) for Offloading and Computation
    while 1

    cvx_begin
        variable x(U, N, M) nonnegative
        variable fUser(U, N) nonnegative
        variable fUAV(U, N, M) nonnegative
        variable threshold
        variable userRate(U, N, M)
        eta = 0;
        expressions pUser(U, N)

        for i = 1:U
            for n = 1:N
                    pUser(i, n) = 0;
            end
        end


        for i = 1:U
            for n = 1:N
                for m = 1:M
                    if Psi(i, n, m) == 1
                        pUser(i, n) = pUser(i, n)+B*N0/g(i, n, m)*(exp(J*x(Psi(i, n, m), n, m))-1);
                    elseif Psi(i, n, m) > 1
                        pUser(i, n) = pUser(i, n)+B*N0/g(i, n, m)*(...
                            exp(J*x(Psi(i, n, m), n, m))-...
                            (...
                            exp(J*x0(Psi(i, n, m)-1, n, m))+J*exp(J*x0(Psi(i, n, m)-1, n, m))*(x(Psi(i, n, m)-1, n, m)-x0(Psi(i, n, m)-1, n, m))...
                            )...
                            );
                    end
                end
            end
        end

        for i = 1:U
            for n = 1:N
                for m = 1:M
                    if Psi(i, n, m) > 1
                        userRate(i, n, m) == x(Psi(i, n, m), n, m)-x(Psi(i, n, m)-1, n, m);
                    elseif Psi(i, n, m) == 1
                        userRate(i, n, m) == x(Psi(i, n, m), n, m);
                    end                
                end
            end
        end

        for m = 1:M
            eta = eta+c2*sum(sum(mu(2:end, :, m).*exp(J*x(:, :, m)))) +...
                c2*sum(mu(1, :, m));
        end
        eta = eta+c1*sum(sum(pow_p(fUser, 3))) +...
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
            for k = 2:U
                for n = 1:N
                    for m = 1:M
                        x(k, n, m) >= x(k-1, n, m);
                    end              
                end
            end

            pUser <= pmax;

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
    break
            
            if abs(eta-eta0) <= 1
                break
            end

            x0 = x;
            eta0 = eta;
    end

%% Iteration begins
while 1    
    %% Optimization of UAV Trajectory Scheduling
    for m = 1:M
        A(:, :, m) = exp(J*x(:, :, m));
    end
    chi = zeros(U, N, M);
    for m= 1:M
        chi(1, :, m) = A(1, :, m)-ones(1, N);
        for k = 2:U
            chi(k, :, m) = A(k, :, m)-A(k-1, :, m);
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
    chi = chi*c4;

    % nu is a U*N*2 matrix.
    % nu(k, n, :) gives coordinates of SMD with k-th smallest channel at time n.
    nu = zeros(U, N, 2);
    for n = 1:N
        for k = 1:U
            nu(k, n, :) = rUser(rank(k, n), 1:2);
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
            for k = 1:U
                for m = 1:M
                Eoffload = Eoffload+chi(k, n, m).*...
                    (...
                    p1*square_pos...
                    (...
                    norm([r1(n, m),r2(n, m),hUAV]-[nu(k, n, 1), nu(k, n, 2), 0])...
                    )...
                    +p2*norm([r1(n, m),r2(n, m),hUAV]-[nu(k, n, 1), nu(k, n, 2), 0])...
                    +p3...
                    );
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

    gRank = zeros(U, N, M);
    rank = zeros(U, N, M);
    Psi = zeros(U, N, M);
    beta = zeros(U+2, N, M);
    mu = zeros(U+1, N, M);

    for m = 1:M
        [gRank(:, :, m), rank(:, :, m), Psi(:, :, m)] = chanRank(g(:, :, m), U, N);

        beta(2:(U+1), :, m) = 1./gRank(:, :, m);

        for k = 1:U+1
            mu(k, :, m) = -beta(k+1, :, m)+beta(k, :, m);
        end
    end
    % rank gives a list of user index. Its value denotes the user index with that ranking, 
    % i.e. rank(k, n) means the index of user with the k-th smallest channel gain at time n.
    % Psi helps track the rankings of users, 
    % i.e. Psi(i, n) denotes the ranking of user i at time n.
    % beta is a (U+2)*N matrix.
    mu = mu/(1e9);

    % Set initial values for iteration.
    x0 = zeros(U, N, M);
    eta0 = 0;

    %% Resource Allocation (RA) for Offloading and Computation
    while 1

    cvx_begin
        variable x(U, N, M) nonnegative
        variable fUser(U, N) nonnegative
        variable fUAV(U, N, M) nonnegative
        variable threshold
        variable userRate(U, N, M)
        eta = 0;
        expressions pUser(U, N)

        for i = 1:U
            for n = 1:N
                    pUser(i, n) = 0;
            end
        end


        for i = 1:U
            for n = 1:N
                for m = 1:M
                    if Psi(i, n, m) == 1
                        pUser(i, n) = pUser(i, n)+B*N0/g(i, n, m)*(exp(J*x(Psi(i, n, m), n, m))-1);
                    elseif Psi(i, n, m) > 1
                        pUser(i, n) = pUser(i, n)+B*N0/g(i, n, m)*(...
                            exp(J*x(Psi(i, n, m), n, m))-...
                            (...
                            exp(J*x0(Psi(i, n, m)-1, n, m))+J*exp(J*x0(Psi(i, n, m)-1, n, m))*(x(Psi(i, n, m)-1, n, m)-x0(Psi(i, n, m)-1, n, m))...
                            )...
                            );
                    end
                end
            end
        end

        for i = 1:U
            for n = 1:N
                for m = 1:M
                    if Psi(i, n, m) > 1
                        userRate(i, n, m) == x(Psi(i, n, m), n, m)-x(Psi(i, n, m)-1, n, m);
                    elseif Psi(i, n, m) == 1
                        userRate(i, n, m) == x(Psi(i, n, m), n, m);
                    end                
                end
            end
        end

        for m = 1:M
            eta = eta+c2*sum(sum(mu(2:end, :, m).*exp(J*x(:, :, m)))) +...
                c2*sum(mu(1, :, m));
        end
        eta = eta+c1*sum(sum(pow_p(fUser, 3))) +...
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
            for k = 2:U
                for n = 1:N
                    for m = 1:M
                        x(k, n, m) >= x(k-1, n, m);
                    end              
                end
            end

            pUser <= pmax;

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

            if abs(eta-eta0) <= 1
                break
            end

            x0 = x;
            eta0 = eta;
    end

    %% Update the weighted system 
    % Update the each term of weighted system utility.
    E_user_offload = 0;
    for m = 1:M
        E_user_offload = E_user_offload+c2*sum(sum(mu(2:end, :, m).*exp(J*x(:, :, m)))) +...
            c2*sum(mu(1, :, m));
    end

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

% Record the offloading power of each user in each time slot.
recordOffloadPower = zeros(U, N);
for i = 1:U
    for n = 1:N
        for m = 1:M
            if Psi(i, n, m) > 1
                recordOffloadPower(i, n) = recordOffloadPower(i, n)+...
                    B*N0/g(i, n, m)*(exp(J*x(Psi(i, n, m), n, m))-exp(J*x(Psi(i, n, m)-1, n, m)));
            elseif Psi(i, n, m) == 1
                recordOffloadPower(i, n) = recordOffloadPower(i, n)+...
                    B*N0/g(i, n, m)*(exp(J*x(Psi(i, n, m), n, m))-1);
            end
        end
    end
end

power_NOMA = recordOffloadPower;
rUAV_NOMA = rUAV;
% utilitySeq_NOMA = utilitySeq(1: numIteration);
fUAV_NOMA = fUAV;
fUser_NOMA = fUser;
utility_NOMA = [E_user_offload, E_user_comp, E_uav_comp, E_uav_move];


end

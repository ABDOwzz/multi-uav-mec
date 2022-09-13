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
% rUser = [lengthArea*rand(U, 1),widthArea*rand(U, 1), zeros(U, 1)];
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
hUAV = 20;                                         % Flying altitude of UAV
rI = [0, 0, hUAV; 0, 0, hUAV];              % Intial position of the UAV (M*3)
rF = [lengthArea, 0, hUAV; 0, widthArea, hUAV];               % Final position of the UAV

rUAV = zeros(N, 3, M);                             % Trajectory vector of the UAV

% rUAV is a N*3 matrix of which each describes the 3-D coordinates of UAV
% at each time slot.
for m = 1:M
    rUAV(N, :, m) = rF(m, :);
end
for n = 1:N-1
    for m = 1:M
        rUAV(n, :, m) = rI(m, :)+n/N*(rF(m, :)-rI(m, :));
    end
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

% Find the channel gain for each user in each time slot
g = zeros(U, N, M);                             % Channel gain matrix
L = zeros(U, N, M);                             % Path loss matrix
% g is a U*N matrix denoting the channel gains for all links throughout the interval.

for i = 1:U
    for n = 1:N
        for m = 1:M
            [g(i, n, m), L(i, n, m)]  = chanGainFSPL(rUAV(n, :, m), rUser(i, :));
        end
    end
end

omegaUser = 0.8;
omegaUAV = 1-omegaUser;

N0 = 1e-17;                     % Power spectral density of AWGN (W/Hz)
B = 4e6;                                        % System bandwidth (Hz)

D = 120*ones(U, 1);                     % Number of bits for each SMDs to finish (Mbit)
kappa = 1e-28;                            % Effective swiched capacitance of CPU

fmaxUAV = 20;                            % Maximum CPU frequency of MEC server (GHz)
fminUAV = 0;                               % Minimum CPU frequency of SMDs (GHz)
CUAV = 1e3;                                % Processing density of MEC server (Hz/bit)

fmaxUser = 1;                          % Maximum CPU frequency of SMDs (GHz)
fminUser = 0;                              % Minimum CPU frequency of SMDs (GHz)
CUser = 1e3;                               % Processing density of SMD (Hz/bit)

%% Resource Allocation (RA) for Offloading and Computation
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
J = log(2)/(B*tau);
J = J*1e6;
c1 = omegaUser*kappa*tau*1e27;
c2 = omegaUser*B*N0*1e9;
c3 = omegaUAV*kappa*tau*1e27;

% Solve the convex optimization problem: Posynomial + Sum of exp

xmin = 0;
cvx_begin
    variable x(U, N, M) nonnegative
    variable fUser(U, N) nonnegative
    variable fUAV(U, N, M) nonnegative
    variable threshold
    variable userRate(U, N, M)
    eta = 0;
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

        % C3: Offloading rate constraint on SMDs
        xmin <= x;
        
        for k = 2:U
            for n = 1:N
                for m = 1:M
                    x(k, n, m) >= x(k-1, n, m);
                end              
            end
        end
        
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

% Record the offloading power of each user in each time slot.
recordPowerOffload = zeros(U, N, M);
for i = 1:U
    for n = 1:N
        for m = 1:M
            if Psi(i, n, m) > 1
                recordPowerOffload(i, n, m) = B*N0/g(Psi(i, n, m))*(exp(J*x(Psi(i, n, m), n, m))-exp(J*x(Psi(i, n, m)-1, n, m)));
            elseif Psi(i, n, m) == 1
                recordPowerOffload(i, n, m) = B*N0/g(Psi(i, n, m))*(exp(J*x(Psi(i), n, m))-1);
            end
        end
    end
end

figure(2)
for i = 1:U
    for m = 1:M
            stem3(i*ones(N, 1), (1:N), recordPowerOffload(i, :, m));
    hold on
    end
end
xlabel('User Index')
ylabel('Time Slot')
zlabel('Offloading Power')


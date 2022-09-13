% clear all
% close all
% lengthArea = 1000;
% widthArea = 1000;
% U = 8;
% rUser = [0.3*lengthArea+0.7*lengthArea*rand(U, 1),0.3*widthArea+0.7*widthArea*rand(U, 1), zeros(U, 1)];

% rUser = [kron(100*(2:2:8).',ones(2,1)), ...
%     repmat([600, 800].', 4, 1),...
%     zeros(U, 1)];

% [power_OMA, rUAV_OMA, utility_OMA, fUAV_OMA, fUser_OMA] = experiment0_OMA(rUser);
[power_NOMA, rUAV_NOMA, utility_NOMA, fUAV_NOMA, fUser_NOMA] = experiment0_NOMA(rUser);

E_user_offload_NOMA =  utility_NOMA(1);
E_user_comp_NOMA = utility_NOMA(2);
E_uav_comp_NOMA = utility_NOMA(3);
E_uav_move_NOMA = utility_NOMA(4);

E_user_offload_OMA =  utility_OMA(1);
E_user_comp_OMA = utility_OMA(2);
E_uav_comp_OMA = utility_OMA(3);
E_uav_move_OMA = utility_OMA(4);

omegaUser = 0.8;                        % Weight factor of users
omegaUAV = 1-omegaUser;         % Weight factor of UAVs


y = [E_user_offload_NOMA/omegaUser, E_user_offload_OMA/omegaUser;...
    E_user_comp_NOMA/omegaUser, E_user_comp_OMA/omegaUser;...
    E_uav_comp_NOMA/omegaUAV, E_uav_comp_OMA/omegaUAV;...
    E_uav_move_NOMA/omegaUAV, E_uav_move_OMA/omegaUAV];

E_total_NOMA = sum(utility_NOMA)
E_total_OMA = sum(utility_OMA)

N = 20;

figure(1)
scatter(rUser(:, 1), rUser(:, 2));
hold on
scatter(rUAV_NOMA(:, 1, 1), rUAV_NOMA(:, 2, 1));
hold on
scatter(rUAV_NOMA(:, 1, 2), rUAV_NOMA(:, 2, 2));
hold on
scatter(rUAV_OMA(:, 1, 1), rUAV_OMA(:, 2, 1));
hold on
scatter(rUAV_OMA(:, 1, 2), rUAV_OMA(:, 2, 2));
legend('Users', 'UAV 1 for NOMA', 'UAV 2 for NOMA', 'UAV 1 for OMA', 'UAV 2 for OMA')

title('Location of SMDs and Trajectories of UAVs')
xlabel('Length/m')
ylabel('Width/m')
axis([0 lengthArea 0 widthArea])

figure(2)
for i = 1:U
	stem3(i*ones(N, 1), (1:N), power_NOMA(i, :));
    hold on
end
xlabel('User Index')
ylabel('Time Slot')
zlabel('Offloading Power')
title('Record of offloading power of each user for NOMA')

figure(3)
for i = 1:U
	stem3(i*ones(N, 1), (1:N), power_OMA(i, :));
    hold on
end
xlabel('User Index')
ylabel('Time Slot')
zlabel('Offloading Power')
title('Record of offloading power of each user for OMA')

figure(4)
bar(y);




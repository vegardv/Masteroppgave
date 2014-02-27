
clc;
clf;
clear;

% MODEL
% Parameters
m = 0.468;
I = diag([4.856e-3,4.856e-3,8.801e-3]);%diag([4.856e-3, 4.856e-3, 8.801e-3]);
g = 9.81;
b = 1.140e-7;
k = 2.980e-6;
l = 0.225;
time = 10;

% CONTROLLER
% Controller attitude
K_p = diag([0.25,0.25,0.25]);
K_i = diag([0,0,0]);
K_d = diag([-0.25,-0.25,-0.25]);

% Controller altitude
k_p = 10;
k_i = 10;
k_d = -5;

% Controller position
K_p_pos = 0.1*eye(2);
K_d_pos = 2*sqrt(K_p_pos(1,1))*eye(2);


% INITIAL VALUES
eta_0 = [0,0,-10,0,0,0]';
nu_0 = [0,0,0,0,0,0]';

% REFERENCE VALUES
psi_ref = 0;
pos_ref = [10,20,-10]';

% SIMULATE
sim('test');


% PLOT RESULTS
% Plot pos
figure(1)
subplot(3,1,1)
plot(UAVPos.Time, UAVPos.Data(1:end, 1), UAVPos.Time, pos_ref.Data(1:end, 1));
ylabel('N [m]');
title('Position of UAV');
legend('Position', 'Reference')


subplot(3,1,2)
plot(UAVPos.Time, UAVPos.Data(1:end, 2), UAVPos.Time, pos_ref.Data(1:end, 2));
ylabel('E [m]');

subplot(3,1,3)
plot(UAVPos.Time, UAVPos.Data(1:end, 3), UAVPos.Time, D_ref.Data(1:end, 1));
xlabel('time [s]');
ylabel('D [m]');

% Plot attitude
figure(2)

subplot(3,1,1)
plot(UAVAttitude.Time, UAVAttitude.Data(1:end, 1), UAVPos.Time, Theta_ref.Data(1:end, 1));
ylabel('\phi [rad]');
title('Attitude of UAV');
legend('Attitude', 'Reference')

subplot(3,1,2)
plot(UAVAttitude.Time, UAVAttitude.Data(1:end, 2), UAVPos.Time, Theta_ref.Data(1:end, 2));
ylabel('\theta [rad]');

subplot(3,1,3)
plot(UAVAttitude.Time, UAVAttitude.Data(1:end, 3), UAVPos.Time, Theta_ref.Data(1:end, 3));
xlabel('time [s]');
ylabel('\psi [rad]');

% Plot ned XY
% figure(3)
% plot(UAVPos.Data(1:end,2), UAVPos.Data(1:end,1));

% Plot u
figure(4)
subplot(2,2,1)
plot(u.Time, u.Data(1:end, 1));

subplot(2,2,2)
plot(UAVAttitude.Time, u.Data(1:end, 2));

subplot(2,2,3)
plot(UAVAttitude.Time, u.Data(1:end, 3));
title('u')

subplot(2,2,4)
plot(UAVAttitude.Time, u.Data(1:end, 4));

% Plot omega
figure(5)
subplot(3,2,1)
plot(u.Time, omega.Data(1:end, 1));

subplot(3,2,2)
plot(UAVAttitude.Time, omega.Data(1:end, 2));

subplot(3,2,3)
plot(UAVAttitude.Time, omega.Data(1:end, 3));
title('omega')

subplot(3,2,4)
plot(UAVAttitude.Time, omega.Data(1:end, 4));

subplot(3,2,5)
plot(UAVAttitude.Time, omega.Data(1:end, 5));
title('omega')

subplot(3,2,6)
plot(UAVAttitude.Time, omega.Data(1:end, 6));

figure(6)
pathPlotterNode(UAVPos.Data(1:end,1), UAVPos.Data(1:end,2),  cornersX.Data(:,:), cornersY.Data(:,:), UAVPos.Time(2), 10, inframe.Data(:), [10,10]');
%pathPlotter(UAVPos.Data(1:end,1), UAVPos.Data(1:end,2), UAVAttitude.Data(1:end,3), UAVPos.Time(2), 80, 0, time, 0,0);
%pathPlotter3d(UAVPos.Data(1:end,1), UAVPos.Data(1:end,2),UAVPos.Data(1:end,3), UAVAttitude.Data(1:end,3), UAVPos.Time(2), 80);
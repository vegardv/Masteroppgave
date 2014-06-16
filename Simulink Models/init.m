
clc;
clf;
close all;
clear;

% MODEL
% Parameters
m = 3;
I = diag([0.1, 0.1, 0.2]);
g = 9.81;
b = 0.1;
k = 3;
l = 0.25;
time = 400;


% CONTROLLER
% Controller attitude
K_p = diag([20,20,0.45]);
K_i = diag([0,0,0]);
K_d = diag([-1,-1,-0.5]);

% Controller altitude
k_p = 100;
k_i = 150;
k_d = -20;

% Controller position
K_p_pos = 0.5*eye(2);
K_d_pos = 2*sqrt(K_p_pos(1,1))*eye(2);
K_i_pos = 0.05 * eye(2);

% INITIAL VALUES
eta_0 = [0,0,0,0,0,0]';
nu_0 = [0,0,0,0,0,0]';

% Node position
nodePosition = [0.5,0.5,0]';

% Disturbance
wind = [0.5,0]';

% REFERENCE VALUES
psi_ref = 0;
pos_ref = [10,10,-10]';

% SIMULATE
sim('test');

% PLOT RESULTS
% Plot pos
figure(1)
subplot(3,1,1)
plot(UAVPos.Time, UAVPos.Data(1:end, 1), 'b', UAVPos.Time, pos_ref.Data(1:end, 1), 'r');
ylabel('N [m]');
legend('Position', 'Reference')
%ylim([-0.2,0.6]);
grid on

subplot(3,1,2)
plot(UAVPos.Time, UAVPos.Data(1:end, 2), 'b', UAVPos.Time, pos_ref.Data(1:end, 2), 'r');
ylabel('E [m]');
grid on
subplot(3,1,3)
plot(UAVPos.Time, UAVPos.Data(1:end, 3), 'b', UAVPos.Time, D_ref.Data(1:end, 1), 'r');
xlabel('time [s]');
ylabel('D [m]');
grid on
%saveas(1, 'C:\Users\vegardvo\Dropbox\Masteroppgave\Report\fig\plots\simulation\positionVarDisturbance2.eps', 'eps2c')


% Plot attitude
figure(2)
subplot(3,1,1)
plot(UAVAttitude.Time, UAVAttitude.Data(1:end, 1), 'b', UAVPos.Time, Theta_ref.Data(1:end, 1), 'r');
ylabel('\phi [rad]');
legend('Attitude', 'Reference')
grid on

subplot(3,1,2)
plot(UAVAttitude.Time, UAVAttitude.Data(1:end, 2), 'b', UAVPos.Time, Theta_ref.Data(1:end, 2), 'r');
ylabel('\theta [rad]');
grid on

subplot(3,1,3)
plot(UAVAttitude.Time, UAVAttitude.Data(1:end, 3),'b', UAVPos.Time, Theta_ref.Data(1:end, 3), 'r');
grid on
xlabel('time [s]');
ylabel('\psi [rad]');
%saveas(2, 'C:\Users\vegardvo\Dropbox\Masteroppgave\Report\fig\plots\simulation\AttNoDisturbance2.eps', 'eps2c')

% Plot u
figure(4)
subplot(2,2,1)
plot(u.Time, u.Data(1:end, 1));
xlabel('time [s]');
ylabel('T [N]')
grid on
subplot(2,2,2)
plot(UAVAttitude.Time, u.Data(1:end, 2));
xlabel('time [s]');
ylabel('\tau_\phi [Nm]');
grid on
subplot(2,2,3)
plot(UAVAttitude.Time, u.Data(1:end, 3));
xlabel('time [s]');
ylabel('\tau_\theta [Nm]');
grid on
subplot(2,2,4)
plot(UAVAttitude.Time, u.Data(1:end, 4));
xlabel('time [s]');
ylabel('\tau_\psi [Nm]');
grid on
%saveas(4, 'C:\Users\vegardvo\Documents\GitHub\Masteroppgave\Report\fig\plots\simulation\uNodeDisturbance.eps', 'eps2c')

% Plot omega
figure(5)
subplot(3,2,1)
plot(u.Time, omega.Data(1:end, 1));
xlabel('time [s]');
ylabel('\omega_1 [rad/s]');
grid on
subplot(3,2,2)
plot(UAVAttitude.Time, omega.Data(1:end, 2));
xlabel('time [s]');
ylabel('\omega_2 [rad/s]');
grid on
subplot(3,2,3)
plot(UAVAttitude.Time, omega.Data(1:end, 3));
xlabel('time [s]');
ylabel('\omega_3 [rad/s]');
grid on
subplot(3,2,4)
plot(UAVAttitude.Time, omega.Data(1:end, 4));
xlabel('time [s]');
ylabel('\omega_4 [rad/s]');
grid on
subplot(3,2,5)
plot(UAVAttitude.Time, omega.Data(1:end, 5));
xlabel('time [s]');
ylabel('\omega_5 [rad/s]');
grid on
subplot(3,2,6)
plot(UAVAttitude.Time, omega.Data(1:end, 6));
xlabel('time [s]');
ylabel('\omega_6 [rad/s]');
grid on
%saveas(5, 'C:\Users\vegardvo\Documents\GitHub\Masteroppgave\Report\fig\plots\simulation\omegaNodeDisturbance.eps', 'eps2c')


figure(6)
pathPlotterNode(UAVPos.Data(1:end,1), UAVPos.Data(1:end,2),  cornersX.Data(:,:), cornersY.Data(:,:), UAVPos.Time(2), 70, inframe.Data(:), nodePos.Data(:, 1:2));
grid on
%saveas(6, 'C:\Users\vegardvo\Dropbox\Masteroppgave\Report\fig\plots\simulation\positionFrameVarDisturbance2.eps', 'eps2c')

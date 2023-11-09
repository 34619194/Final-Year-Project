%Author: Oskar Lund,
%We Description
%11.10. 17
% This script simulates the overall behaviour of the
%butterfly robot in open loop, in addition to new,
%potential trajectories

%Define the object BR of class Butterfly Robot
BR = ButterflyRobot;

%We Choose which program to run:ABG_EOM
program = 1;
interval = 2.5 ;
if (program == 0)
%Simulate alpha, beta ,

t0_  = 0; %Initail time 
tEnd_ = 3; %final time
x0_ = [0; 4.3];
[t_,x_] = BR.simABG(t0_ , tEnd_ , x0_);
plot(t_,x_)

elseif (program ==1)
%Simulate the full butterfly system

t0 = 0;
tEnd =interval;
x0 = [0; 0; 0.086; 4.3] ;
[t ,x ,Fn , Fs , u] = BR.simBR(t0, tEnd, x0);
% 1- angle of frame, 2- angle of ball, 3- deriviatave angel of frame, 4 -
% deriviative of angle of ball

end 
% Create a single figure with subplots
frame_angle = x(:,1)./2;
ball_angle = x(:,2)./2;
dframe_angle = x(:,3)./4;
dball_angle = x(:,4)./4;
acceleration  = derivative(t,dframe_angle);
torque_all = acceleration.*BR.Jf;
theta = max(frame_angle);
Dtheta = max(dframe_angle)
DDtheta = max(acceleration);
torque = DDtheta*BR.Jf

figure;

% Plot 1: Angle
subplot(2, 2, 1);
plot(t, frame_angle);
hold on;
plot(t, ball_angle);
xlim([0 interval]);
ylim([0 2 * pi()]);
title("Simulated ${\theta}$ and ${\varphi}$ values", 'interpreter', 'latex');
legend("$\theta (t)$", "$\varphi (t)$", 'interpreter', 'latex');
yticks([0, 1 / 2 * pi, pi, 3 / 2 * pi, 2 * pi]);
yticklabels({'0', '1/2\pi', '\pi', '3/2\pi', '2\pi'});
xlabel("Time [s]");
ylabel("Angle [rad]");
grid on;
hold off;

% Plot 2: Angular Velocity
subplot(2, 2, 2);
plot(t, dframe_angle);
hold on;
plot(t, dball_angle);
xlim([0 interval]);
grid on;
legend("$\dot\theta (t)$", "$\dot\varphi (t)$", 'interpreter', 'latex');
title("Simulated ${\dot\theta}$ and ${\dot\varphi}$ values", 'interpreter', 'latex');
yticks([-2 * pi, -3 / 2 * pi, -pi, -1 / 2 * pi, 0, (1 / 2) * pi, pi, 3 / 2 * pi, 2 * pi, 5 / 2 * pi]);
yticklabels({'-2\pi', '-3/2\pi', '-\pi', '-1/2\pi', '0', '1/2\pi', '\pi', '3/2\pi', '2\pi', '5/2\pi'});
xlabel("Time [s]");
ylabel("Angular Velocity [rad/s]");
hold off;

% Plot 3: Acceleration
subplot(2, 2, 3);
plot(t, acceleration);
grid on;
title("Acceleration vs time", 'interpreter', 'latex');
xlabel("Time [s]");
ylabel("Acceleration [rad/s^2]");
xlim([0 interval]);

% Plot 4: Torque
subplot(2, 2, 4);
plot(t, torque_all);
grid on;
title("Torque vs time", 'interpreter', 'latex');
xlabel("Time [s]");
ylabel("Torque [Nm]");
xlim([0 interval]);



 %clear all; close all;
 close all; clear; clc; close all;


%% Double Pendulum System Parameters
gv = 9.81; % Gravity Constant
m1 = 1;    % Mass of bob 1 
m2 = 1;    % Mass of bob 2
L1 = 1;    % Length of Rod 1
L2 = 1/4;  % Length of Rod 2 (Cases: [1/4 1/2] m)

%% Initial Conditions
ICs = [pi/3 0 pi*2/3 0];

%% Time Step 
t_intv = [0 5];
t_step = 500;

%% Integrators

[t_rk,y_rk]   = pend(t_intv,ICs,t_step,m1,m2,L1,L2,gv,'rk');
[t_trp,y_trp] = pend(t_intv,ICs,t_step,m1,m2,L1,L2,gv,'trp');
[t_eul,y_eul] = pend(t_intv,ICs,t_step,m1,m2,L1,L2,gv,'eul');

%% ODE 45 

[T,Y] = ode45(@YDOT,t_intv,ICs);

%% Plot

figure
hold on

% Plot data for different methods
plot(t_eul, theta_eul, ':', 'color', "#0072BD", 'LineWidth', 2)
plot(t_trp, theta_trp, '-.', 'color', "#D95319", 'LineWidth', 2)
plot(t_rk, theta_rk, '--', 'color', "#EDB120", 'LineWidth', 2)
plot(t_ode45, theta_ode45, 'k-', 'LineWidth', 2)

xlabel ('Time [s]')
ylabel ('\theta_1 [rad]')
legend('Euler','Trapezoid','Runge-Kutta','ODE45','Location','best')
ylim([-5 5])
grid on
box on 
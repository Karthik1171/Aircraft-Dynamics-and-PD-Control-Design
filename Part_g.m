clc; clear; close all;

m = 0.1;    
Mc = 1;     
L = 0.5;    
g = 9.81;   

A = [0, 1; ((Mc + m) * g) / (L * Mc), 0];
B = [0; -1 / (L * Mc)];
C = [1, 0];
D = 0;

sys = ss(A, B, C, D);

% PD controller design
zeta = 0.5913;  
wn = 4.4 / (zeta * 1);  

Kp = -(wn^2 + 21.582)/2;
Kd = -zeta * wn;

% Closed-loop system
A_cl = [0, 1; ((Mc + m) * g + Kp) / (L * Mc), Kd / (L * Mc)];
sys_cl = ss(A_cl, B, C, D);

% Time span
tspan = linspace(0, 10, 500);
dt = tspan(2) - tspan(1);

% Initial conditions
theta0 = pi/18;  
theta_dot0 = 0;
y0 = [theta0; theta_dot0];

% Simulate linearized closed-loop system
[y_lin, t, x_lin] = initial(sys_cl, y0, tspan);

% Nonlinear system simulation using RK4
x_nonlin = zeros(length(tspan), 2);
x_nonlin(1, :) = y0;

for j = 1:length(tspan) - 1
    t_j = tspan(j);
    x_j = x_nonlin(j, :)';
    
    k1 = nonlinear_pendulum(t_j, x_j, m, Mc, L, g, Kp, Kd);
    k2 = nonlinear_pendulum(t_j + dt/2, x_j + dt/2 * k1, m, Mc, L, g, Kp, Kd);
    k3 = nonlinear_pendulum(t_j + dt/2, x_j + dt/2 * k2, m, Mc, L, g, Kp, Kd);
    k4 = nonlinear_pendulum(t_j + dt, x_j + dt * k3, m, Mc, L, g, Kp, Kd);
    
    x_nonlin(j+1, :) = (x_j + (dt/6) * (k1 + 2*k2 + 2*k3 + k4))';
end

% Plot results
figure;
subplot(2,1,1);
plot(t, x_lin(:, 1) , 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\theta (rad)');
title('PD Controlled Linearized Response');
grid on;

subplot(2,1,2);
plot(t, x_nonlin(:, 1) , 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\theta (rad)');
title('Original Dynamical System Response ');
grid on;

% Nonlinear system equations
function dxdt = nonlinear_pendulum(~, x, m, Mc, L, g, Kp, Kd)
    theta = x(1);
    theta_dot = x(2);
    
    % PD control law
    u = -Kp * (theta) - Kd * theta_dot ;

    % Nonlinear dynamics
    F = (sin(theta) / (L * (Mc + m * sin(theta)^2))) * (-m * L * theta_dot^2 * cos(theta) + (Mc + m) * g);
    G = (-cos(theta)) / (L * (Mc + m * sin(theta)^2));
    
    theta_ddot = F + G * u;
    
    dxdt = [theta_dot; theta_ddot];
end

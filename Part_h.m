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
zeta = 0.517;  
wn = 4.4/ (zeta * 0.5);  

Kp = (-wn^2 - 21.582)/2;
Kd = -zeta * wn ;

% Closed-loop linearized system
A_cl = [0, 1; ((Mc + m) * g + Kp) / (L * Mc), Kd / (L * Mc)];
sys_cl = ss(A_cl, B, C, D);

tspan = linspace(0, 10, 500);
dt = tspan(2) - tspan(1); 

theta0 = pi/18;  
theta_dot0 = 0;
y0 = [theta0; theta_dot0];

[y_lin, t, x_lin] = initial(sys_cl, y0, tspan);

x_nonlin = zeros(length(tspan), 2);
x_nonlin(1, :) = y0;

for i = 1:length(tspan) - 1
    t_i = tspan(i);
    x_i = x_nonlin(i, :)';
    
    k1 = nonlinear_pendulum(t_i, x_i, m, Mc, L, g, Kp, Kd);
    k2 = nonlinear_pendulum(t_i + dt/2, x_i + dt/2 * k1, m, Mc, L, g, Kp, Kd);
    k3 = nonlinear_pendulum(t_i + dt/2, x_i + dt/2 * k2, m, Mc, L, g, Kp, Kd);
    k4 = nonlinear_pendulum(t_i + dt, x_i + dt * k3, m, Mc, L, g, Kp, Kd);
    
    x_nonlin(i+1, :) = (x_i + (dt/6) * (k1 + 2*k2 + 2*k3 + k4))';
end

% Plot results
figure;
subplot(2,1,1);
plot(t, x_lin(:, 1), 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\theta (rad)');
title('PD Controlled Linearized Response');
grid on;

subplot(2,1,2);
plot(t, x_nonlin(:, 1), 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\theta (rad)');
title('Original Dynamical System Response ');
grid on;

% Nonlinear system equations
function dxdt = nonlinear_pendulum(t, x, m, Mc, L, g, Kp, Kd)
    theta = x(1);
    theta_dot = x(2);
    F = (sin(theta) / (L * (Mc + m * sin(theta)^2))) * (-m * L * theta_dot^2 * cos(theta) + (Mc + m) * g);
    G = (-cos(theta)) / (L * (Mc + m * sin(theta)^2));
    
    % PD control law
    u = -Kp * (theta) - Kd * theta_dot;
    
    dxdt = zeros(2, 1);
    dxdt(1) = theta_dot;
    dxdt(2) = F + G * u;
end

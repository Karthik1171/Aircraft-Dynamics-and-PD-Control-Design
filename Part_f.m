clc; clear; close all;

m = 0.1;    
Mc = 1;     
L = 0.5;    
g = 9.81;   

A_matrices = {
    [0, 1; ((Mc + m) * g) / (L * Mc), 0],  
    [0, 1; -((Mc + m) * g) / (L * Mc), 0]  
};
B_matrices = {
    [0; -1 / (L * Mc)],  
    [0; 1 / (L * Mc)]    
};

C = [1, 0];  
D = 0;  

tspan = linspace(0, 20, 500);  
dt = tspan(2) - tspan(1);  

theta_eq = [0, pi];  

figure;
tiledlayout(2, 2); 

for i = 1:2
    A = A_matrices{i};
    B = B_matrices{i};
    sys = ss(A, B, C, D);

    theta0 = theta_eq(i) + pi/18; 
    theta_dot0 = 0;
    y0 = [theta0; theta_dot0];

    u = zeros(size(tspan));  
    [y_lin, t, x_lin] = lsim(sys, u, tspan, y0);

    x_nonlin = zeros(length(tspan), 2);
    x_nonlin(1, :) = y0;

    for j = 1:length(tspan) - 1
        t_j = tspan(j);
        x_j = x_nonlin(j, :)';
        
        k1 = pendulumDynamics(t_j, x_j, m, Mc, L, g, 0);
        k2 = pendulumDynamics(t_j + dt/2, x_j + dt/2 * k1, m, Mc, L, g, 0);
        k3 = pendulumDynamics(t_j + dt/2, x_j + dt/2 * k2, m, Mc, L, g, 0);
        k4 = pendulumDynamics(t_j + dt, x_j + dt * k3, m, Mc, L, g, 0);
        
        x_nonlin(j+1, :) = (x_j + (dt/6) * (k1 + 2*k2 + 2*k3 + k4))';
    end
    
    % (nonlinear)
    nexttile;
    plot(t, x_nonlin(:, 1), 'k', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('\theta (rad)');
    title(['Nonlinear response for \theta_{eq} = ' num2str(theta_eq(i))]);
    grid on;

    %  (linear)
    nexttile;
    plot(t, x_lin(:, 1), 'r', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('\theta (rad)');
    title(['Linear response for \theta_{eq} = ' num2str(theta_eq(i))]);
    grid on;
end

function dydt = pendulumDynamics(~, y, m, Mc, L, g, u)
    theta = y(1);
    theta_dot = y(2);

    F = (sin(theta) / (L * (Mc + m * sin(theta)^2))) * (-m * L * theta_dot^2 * cos(theta) + (Mc + m) * g);
    
    G = (-cos(theta)) / (L * (Mc + m * sin(theta)^2));
    
    theta_ddot = F + G * u;

    dydt = [theta_dot; theta_ddot];
end
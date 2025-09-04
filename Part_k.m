clc; clear; close all;

Kp_values = [-10.79, -11.5, -10];  
Kd = 0;  

t = 0:0.01:20;

theta0 = pi/18;  
x0 = [theta0; 0];  

for idx = 1:length(Kp_values)
    Kp = Kp_values(idx);

    f = @(t, x) [
        x(2);  
        21.582 * x(1) - 2 * (-(Kp * x(1)) - Kd * x(2) + 1)  
    ];
    
    x = zeros(2, length(t));
    x(:, 1) = x0;

    for i = 1:length(t)-1
        k1 = f(t(i), x(:, i));
        k2 = f(t(i) + 0.01/2, x(:, i) + 0.01/2 * k1);
        k3 = f(t(i) + 0.01/2, x(:, i) + 0.01/2 * k2);
        k4 = f(t(i) + 0.01, x(:, i) + 0.01 * k3);

        x(:, i+1) = x(:, i) + 0.01/6 * (k1 + 2*k2 + 2*k3 + k4);
    end

    theta = x(1, :);

    subplot(3, 1, idx);
    plot(t, theta, 'r', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (sec)');
    ylabel('\theta (rad)');
    title(['Response for K_p = ', num2str(Kp)]);
end


for idx = 1:length(Kp_values)
    Kp = Kp_values(idx);
    
    num = [2];  
    den = [1, 2 * Kd, 2 * Kp];  
    
    sys = tf(num, den);

    figure;
    pzmap(sys);
    grid on;
    title(['Zeros & Poles Map for K_p = ', num2str(Kp)]);
    xlabel('Real Axis');
    ylabel('Imaginary Axis');
end

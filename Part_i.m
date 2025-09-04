

m = 0.1;    
Mc = 1;     
L = 0.5;    
g = 9.81;  

Kp = -155.65;  
Kd = -8.8;      


f = @(t, x) [
    x(2);  
    21.582 * x(1) - 2 * (-(Kp * x(1)) - Kd * x(2) + 1)  
];

% Time vector
t_start = 0;
t_end = 20;
dt = 0.01; 
t = t_start:dt:t_end;

% Initial conditions 
theta0_rad = pi/18;  
x0 = [theta0_rad; 0];  

x = zeros(2, length(t));
x(:, 1) = x0;

for i = 1:length(t)-1
    k1 = f(t(i), x(:, i));
    k2 = f(t(i) + dt/2, x(:, i) + dt/2 * k1);
    k3 = f(t(i) + dt/2, x(:, i) + dt/2 * k2);
    k4 = f(t(i) + dt, x(:, i) + dt * k3);
    
    x(:, i+1) = x(:, i) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

theta_deg = x(1, :)*(180/pi);

%  steady-state error
e_steady_state = -2 / (2 * Kp + 21.582);
e_steady_state_deg = e_steady_state* (180/pi);
fprintf(['Steady-State Error: ', num2str(e_steady_state_deg), ' degrees']);

figure;
plot(t, theta_deg, 'r', 'LineWidth', 1.5);
yline(e_steady_state_deg, 'k--', 'LineWidth', 1.5);
grid on;
xlabel('Time (sec)');
ylabel('\theta (deg)');
title('Response for Unit Step Disturbance');
legend('Simulated Response', 'Steady-State Error');
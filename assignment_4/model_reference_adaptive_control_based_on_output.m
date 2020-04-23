clear;

% observer roots
Ac = -1;
Bc = 1;
% use Matlab to derive a state space model from the transfer funciton
num = [0,1,1];
denom = [1, -3, 2];
% we note that Dp is 0 for our plant
[Ap, Bp, hp, Dp] = tf2ss(num, denom);
% we get this from the transfer function
kp = 1;
%initial conditions of the model and the plant
ym0 = 0;
xp0 = [0;0];
%initial values of the adaptive parameters
k0 = 0;
theta00 = 0;
theta10 = 0;
theta20 = 0;
%create time span and input command
t_sampling = 0.1;
time = 0:t_sampling:1200;
% r = cos(time) + 10*cos(5*time);
r = randn(length(time), 1);
% call the simulate system function with all necessary parameters to run
% the simulation
zs = simulate_system(ym0, xp0, r, Ap, Bp, hp, Ac, Bc, kp, time,  k0, theta00, theta10, theta20);
% retreive the model output and plant state values over time
yms = zs(:, 1);
xps = zs(:, 2:3);
% retreive the model output and plant state values over time
params = zs(:, 6:9);
% based on plant states calculate the plant output over time
yps = (hp*xps')';
%plot the results
fig = figure;
subplot(2,1,1);
%plot tracking error over time
plot(time, yms- yps);
title("Tracking error over time");
subplot(2,1,2);
%plot parameters over time
plot(time, params);
legend(["k", "\theta_0", "\theta_1", "\theta_2"]);
title("Parameter values over time");

function [zs] = simulate_system(ym0, xp0, r, Ap, Bp, hp, Ac, Bc, kp, time, k0, theta00, theta10, theta20)
    % find the sampling time
    step_size = time(2) - time(1);
    % define the initial conditions for xc1 and xc2
    % equations
    xc10 = 0;
    xc20 = 0;
    % define the collection of initial conditions for the system of
    % dfferential equations
    z0 = [ym0; xp0; xc10; xc20; k0; theta00; theta10; theta20];
    % preallocate space to keep the results
    zs = zeros(length(time), 9);
    for index = 1:(length(time))
        % save the previous results
        zs(index,:) = z0;
        % create the small time span to solve the differential equations
        % for this small period of time
        t_start = time(index);
        r0 = r(index);
        % this is a workaround. The ode45 defines integration steps automatically, and
        % if I specify only t_start and t_stop, it would return an array of
        % about 50 mid-points. With this way it will return only for 3 time
        % points: t_start, (t_start+t_end)/2, t_end)
        tspan = [t_start, t_start + step_size/2, t_start + step_size]; 
        % solve the differential equations for current time span and update
        % the initial conditions to be used for the following time span
        [~, z] = ode45(@(time, z) syst(time, z, r0, Ap, Bp, Ac, Bc, hp, kp), tspan, z0);
        z0 = z(3,:);
    end
end

% this function immitates the system
function dzdt = syst(t, z, r, Ap, Bp, Ac, Bc, hp, kp)
    if false
        if t > 200
            Ap = [2, -1; 1, 0];
        end
    end
    % retrived initial values from z
    ym = z(1);
    xp = z(2:3);
    xc1 = z(4);
    xc2 = z(5);
    k = z(6);
    theta0 = z(7);
    theta1 = z(8);
    theta2 = z(9);
    % as the Dp is zero in our case, we do not need to case about Dp*u term
    % for yp
    % calculate the output of the plant and the control input of the plant
    yp = hp*xp;
    u = k*r + theta0 * yp + theta1*xc1 + theta2*xc2;
    % differential equation for the model and plant simulation
    % as the transfer function of the model was fairly simple, I used the
    % differential equation directly instead of using states first and then
    % calculating the output
    dymdt = -ym + r;
    dxpdt = Ap*xp + Bp*u;
    % differential equations for the xc1 and xc2
    dxc1dt = Ac*xc1 + Bc*u;
    dxc2dt = Ac*xc2 + Bc*yp;
    % calculate the current error value
    e1 = (yp - ym);
    % differential equations for parameter updates
    dkdt =      -sign(kp)*e1*r;
    dtheta0dt = -sign(kp)*e1*yp;
    dtheta1dt = -sign(kp)*e1*xc1;
    dtheta2dt = -sign(kp)*e1*xc2;
    % collect all into one array to return
    dzdt = [dymdt; dxpdt; dxc1dt; dxc2dt; dkdt; dtheta0dt; dtheta1dt; dtheta2dt];

end
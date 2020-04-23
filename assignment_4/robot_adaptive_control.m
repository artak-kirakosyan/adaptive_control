clear
clear global
%define the mass and the length of the robot joints
m1 = 1;
m2 = 2;
l = 1;
% define the model through which we pass the commanded signal
Am = [0, 1; -4, -4];
Bm = [0; 4];
% define the error gains for velocity and position
Kp = 4*eye(2);
Kv = 4*eye(2);
%matrix to be used for Lyapunov calculation
global p Q;
Q = eye(4);
A = [zeros(2,2), eye(2); -Kp, -Kv];
p = sylvester(A', A, -Q);
psi = 4*eye(2);
%%
%initial conditions for desired and actual angles
q1d_0 = [0;0];
q2d_0 = [0;0];
q1_0 = [0;0];
q2_0 = [0;0];
% initial conditions for parameters
P0 = [1.5; 3];
% initial condition used in the ode
z0 = [q1d_0; q2d_0; q1_0; q2_0; P0];

% create time interval and command signals
t_sampling = 0.05;
t = 0:t_sampling:40;
q1_com = square(t/3);
q2_com = square(t/3);

zs = simulate_system(z0, Am, Bm, t, q1_com, q2_com, m1, m2, l, Kp, Kv, psi);

% plot the results
fig = plot_the_results(t, zs);

global V V_dot;
figure
subplot(2,1,1);
plot(1:length(V), V)
title("Lyapunov function(top) and its derivative(bottom)");
subplot(2, 1, 2);
plot(1:length(V_dot), V_dot);

% check if Lyapunov function is valid or not: if it has been increased for
% any time interval, its an invalid one
% if detected invalidity: print the index after which it increased
for i = 1:length(V)-1
    if V(i+1) - V(i) > 0
        disp("********************Lyapunov function is invalid: " + i);
        break;
    end
end

% check if Lyapunov function is valid or not: if it has been increased for
% any time interval, its an invalid one
% if detected invalidity: print the index after which it increased
for i = 1:length(V_dot)
    if V_dot(i) > 0
        disp("********************Lyapunov function derivative is invalid: " + i);
        break;
    end
end


function [zs] = simulate_system(z0, Am, Bm, t, q1_com, q2_com, m1, m2, l, Kp, Kv, psi)
    % find the sampling time
    step_size = t(2) - t(1);
    zs = zeros(length(t), length(z0));
    for index = 1:(length(t))
        % save the previous results
        zs(index,:) = z0;
        % create the small time span to solve the differential equations
        % for this small period of time
        t_start = t(index);
        q1_c = q1_com(index);
        q2_c = q2_com(index);
        % this is a workaround. The ode45 defines integration steps automatically, and
        % if I specify only t_start and t_stop, it would return an array of
        % about 50 mid-points. With this way it will return only for 3 time
        % points: t_start, (t_start+t_end)/2, t_end)
        tspan = [t_start, t_start + step_size/2, t_start + step_size]; 
        % solve the differential equations for current time span and update
        % the initial conditions to be used for the following time span
        [~, z] = ode45(@(t, z) syst(t, z, Am, Bm, q1_c, q2_c, m1, m2, l, Kp, Kv, psi), tspan, z0);
        z0 = z(3,:);
    end
end


% this function immitates the system
function dzdt = syst(t, z, Am, Bm, q1_c, q2_c, m1, m2, l, Kp, Kv, psi)
    if false
        if t > 70
            m2 = 1;
        end
        if t > 140
            m1 = 1.4;
        end
    end
    % first joint's desired position and velocity
    q1d = z(1);
    q1d_dot = z(2);
    z1 = [q1d; q1d_dot];
    % second joint's desired position and velocity
    q2d = z(3);
    q2d_dot = z(4);
    z2 = [q2d; q2d_dot];
    
    % actual join positions and velocities
    q1 = z(5);
    q2 = z(6);
    q1_dot = z(7);
    q2_dot = z(8);
    % rearrange them to use later
    q = [q1;q2];
    x = [q1_dot;q2_dot];
    
    % parameter values
    p1 = z(9);
    p2 = z(10);
    P = [p1; p2];
    % calculate the position and velocity error
    E = [q1d - q1;  q2d - q2];
    E_dot = [q1d_dot - q1_dot; q2d_dot - q2_dot];
    
    % filter the error
    E1 = E_dot + psi * E;
    
    % pass the commanded signal through the model
    dz1dt = Am * z1 + Bm * q1_c;
    dz2dt = Am * z2 + Bm * q2_c;
    
    % separate the joint desired accelerations to be used later
    q1d_ddot = dz1dt(2);
    q2d_ddot = dz2dt(2);
    
    q_star_ddot = [q1d_ddot; q2d_ddot] + Kv * E_dot + Kp*E;
    % determine the M and C(Q in the paper)
    M = [m1*l^2 + 2*m2*l^2 + 2*m2*l^2*cos(q2), m2*l^2 + m2*l^2*cos(q2); m2*l^2 + m2*l^2*cos(q2), m2*l^2];
    C = [-2*m2*l^2*q1_dot*q2_dot*sin(q2) - m2*l^2*q2_dot^2 *sin(q2); m2*l^2*q1_dot^2*sin(q2)];
    
    % determine the M_hat, C_hat and W(Y in paper)
    M_hat = [p1 + 2*p2*(1+cos(q2)), p2*(1+cos(q2)); p2*(1+cos(q2)), p2];
    C_hat = [p2*(-2*q1_dot*q2_dot*sin(q2) - q2_dot^2*sin(q2)); p2*q1_dot^2*sin(q2)];
    
    W = [q_star_ddot(1), 2*q_star_ddot(1) + q_star_ddot(2) + 2*q_star_ddot(1)*cos(q2) + q_star_ddot(2)*cos(q2) - 2*q1_dot*q2_dot*sin(q2) - q2_dot^2*sin(q2); 0, q_star_ddot(1) + q_star_ddot(1)*cos(q2) + q_star_ddot(2) + q1_dot^2*sin(q2)];
    % find the control torque that will be applied to the robot
    % those 2 formulas for T are equivalent
    T = W * P;
    %T = M_hat * q_star_ddot + C_hat;
    
    %  integrate the actual acceleration of the joints(q_ddot)
    % to get the actual velocity and position.   
    % x is the q_dot, then the dxdt will be the q_ddot, which is given as
    % the formula below
    dqdt = x;
    dxdt = M\(T - C);
    
    % adaptation law for parameters
    dPdt = W' * (M_hat\E1);
    % collect all together
    dzdt = [dz1dt; dz2dt; dqdt; dxdt; dPdt];
    
    phi = [1;2] - P;
    global V p V_dot Q;
    V = [V, [E; E_dot]'*p*[E; E_dot] + phi'*phi];
    V_dot = [V_dot, [E; E_dot]'*(-Q)*[E; E_dot] +  2*phi'*W'*inv(M_hat)*E1 - 2*phi'*dPdt];
end

function fig = plot_the_results(t, zs)
    q1d = zs(1:end, 1);
    dq1ddt = zs(1:end, 2);
    q2d = zs(1:end, 3);
    dq2ddt = zs(1:end, 4);

    q1 = zs(1:end, 5);
    q2 = zs(1:end, 6);
    dq1dt = zs(1:end, 7);
    dq2dt = zs(1:end, 8);

    p1 = zs(1:end, 9);
    p2 = zs(1:end, 10);

    fig = figure;
    subplot(3,2,1);
    plot(t, [q1d, q1, q1d-q1]);
    legend(["Desired", "Actual", "Error"]);
    title("First joint desired and actual position");
    subplot(3,2,2);
    plot(t, [dq1ddt, dq1dt, dq1ddt-dq1dt]);
    legend(["Desired", "Actual", "Error"]);
    title("First joint desired and actual velocity");

    subplot(3,2,3);
    plot(t, [q2d, q2, q2d-q2]);
    legend(["Desired", "Actual", "Error"]);
    title("Second joint desired and actual position");
    subplot(3,2,4);
    plot(t, [dq2ddt, dq2dt, dq2ddt-dq2dt]);
    legend(["Desired", "Actual", "Error"]);
    title("Second joint desired and actual velocity");

    subplot(3,2,[5,6]);
    plot(t, [p1, p2]);
    legend(["p1", "p2"]);
    title("Parameter evolution over time");
end
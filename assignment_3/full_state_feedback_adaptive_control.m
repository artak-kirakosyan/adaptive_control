clear;

%define plant's and model's ssytem and input matrices
Ap = [0, 1; -16, -1];
Bp = [0; 1];

Am = [0, 1; -25, -10];
Q = 1;
Bm = Q*Bp;

% solve the equation derived in class to find the p matrix
% Instead of using identity matrix, I used a much bigger one to speed up
% the convergence process
syms p11 p12 p21 p22;
p = [p11, p12; p21, p22];
q = solve(Am'* p + p*Am == -10000*eye(2), {p11, p12, p21, p22});
p = [double(q.p11), double(q.p12); double(q.p21), double(q.p22)];

% solve the equation defined in class to find the ideal values of
% parameters
syms k1 k2;
K_star = [k1, k2];
q = solve(Ap + Bp*K_star == Am, {k1, k2});
K_star= [double(q.k1), double(q.k2)];

% define initial states of the plant, model and the parameters
K0 = [0, 0];
xm0 = [0; 0];
xp0 = [0; 0];

%combine all together to use in ODE solver
z0 = [xm0; xp0; K0'];

%create time span and input command
t_sampling = 10;
time = 0:t_sampling:400;
ucom = randn(length(time), 1);
%ucom = ones(length(time), 1);

% pass everything to the function to simulate the system
[zs, V] = simulate_system(xm0, xp0, K0, K_star, ucom, Am, Bm, Ap, Bp, Q, p, time);

% check if Lyapunov function is valid or not: if it has been increased for
% any time interval, its an invalid one
% if detected invalidity: print the index after which it increased
for i = 1:length(V)-1
    if V(i+1) - V(i) >= 0
        disp("********************Lyapunov function is invalid: " + i);
        break;
    end
end

%%% WARNING: Lyapunv function value is not correct if the plant is dynamic.
%%% Its because I use the old values of the desired values of parameters to
%%% calculate the parameter error, hence the Lyapunov function will be
%%% incorrect

% plot the value of the Lyapunov function
figure;
plot(time, V);
title("Lyapunov function value over time");

% plot the state errors and parameter errors
fig = figure;
subplot(3,1,1);
plot(time, [zs(:,3) - zs(:,1)]);

title("First state error");
legend(["Xp1 - Xm1"]);

subplot(3,1,2);
plot(time, [zs(:,4) - zs(:,2)]);
title("Second state error");
legend(["Xp2 - Xm2"]);

subplot(3,1,3);
plot(time, [zs(:,5), zs(:,6)]);
title("Parameters over time");
legend(["K1", "K2"]);


function [zs, V] = simulate_system(xm0, xp0, K0, K_star, ucom, Am, Bm, Ap, Bp, Q, p, time)
    % find the sampling time
    step_size = time(2) - time(1);
    
    % define the initial conditions for the system of differential
    % equations
    
    z0 = [xm0; xp0; K0'];
    % preallocate space to keep the results
    zs = zeros(length(time), 6);
    V = zeros(length(time), 1);
    for index = 1:(length(time))
        % save the previous results
        zs(index,:) = z0;
        
        % create the small time span to solve the differential equations
        % for this small period of time
        t_start = time(index);
        ucom0 = ucom(index);
        
        % this is a workaround. The ode45 defines integration steps automatically, and
        % if I specify only t_start and t_stop, it would return an array of
        % about 50 mid-points. With this way it will return only for 3 time
        % points: t_start, (t_start+t_end)/2, t_end)
        tspan = [t_start, t_start + step_size/2, t_start + step_size];
        
        % solve the differential equations for current time span and update
        % the initial conditions to be used for the following time span
        [~, curr_z] = ode45(@(t, z) syst(t, z, ucom0, Am, Bm, Ap, Bp, Q, p), tspan, z0);
        z0 = curr_z(3,:);
        
        % take current results, calculate the Lyapunov function and save it
        % to return
        xm = z0(1:2)';
        xp = z0(3:4)';
        e = xp - xm;
        K = z0(5:6);
        phi = K - K_star;

        V(index) = e'*p*e + trace(phi'*phi);
    end
end

% this function is being used by the ode solver
function dzdt = syst(t, z, ucom, Am, Bm, Ap, Bp, Q, p)

    % this part of the code will change the parameters of they plant
    % during the simulation
    % change the condition from false to true to enable it
    if false
        if t > 300
            Ap = [0, 1; -25, -10];
        end
    end
    % take initial values
    xm = z(1:2);
    xp = z(3:4);
    K = z(5:6)';
    
    % calculate the control input to the plant and the state error
    u = K * xp + Q * ucom;
    e = xp - xm;
    
    % differential equations for the model, plant and parameter update
    dxmdt = Am * xm + Bm * ucom;
    dxpdt = Ap * xp + Bp * u;
    dKdt = -Bp'*p*e*xp';
    
    % collect all to return
    dzdt = [dxmdt; dxpdt; dKdt'];
end
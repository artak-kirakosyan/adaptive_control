clear;

% initialize the parameters to be used
am = 4;
bm = 4;
a = 2;
b = 1;

% create the desired and initial values for estimates
td = bm/b;
sd = (am - a)/b;
s0 = sd/2;
t0 = td/2;

% initialize gamma and initial values of the model and plant
gamma = 5;
ym0 = 0;
yp0 = 0;

%create time span and input command
t_sampling = 1;
time = 0:t_sampling:600;

ucom = randn(length(time), 1);
% call the function to simulate the MRAC system

zs = simulate_system(ym0, yp0, ucom, am, bm, a, b, s0, t0, gamma, time);
ym = zs(:,1);
yp = zs(:,2);
s0 = zs(:,3);
t0= zs(:,4);
es = yp - ym;
phi = s0 - sd;
psi = t0 - td;

% calculate energy function value
V = 1/2*(es.^2 + 1/gamma*phi.^2 + 1/gamma*psi.^2);

%%% WARNING: Lyapunv function value is not correct if the plant is dynamic.
%%% Its because I use the old values of the desired values of parameters to
%%% calculate the parameter error, hence the Lyapunov function will be
%%% incorrect

figure;
plot(time, V);
title("Lyapunov function value over time");

fig = figure;
% plot the results
subplot(3,1,1);
plot(time, [ym, yp]);
title("Model response vs plant response. Gamma: " + gamma);
legend(["ym", "yp"]);
subplot(3,1,2);
plot(time, [t0, s0]);
legend(["t0", "s0"]);
title("Evolution of the estimates over time");
subplot(3,1,3);
plot(time, es);
title("Error over time");


function [zs] = simulate_system(ym0, yp0, ucom, am, bm, a, b, s0, t0, gamma, time)
    % find the sampling time
    step_size = time(2) - time(1);
    % define the initial conditions for the system of differential
    % equations
    z0 = [ym0; yp0; s0; t0];
    % preallocate space to keep the results
    zs = zeros(length(time), 4);
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
        [~, z] = ode45(@(time, z) syst(time, z, ucom0, a, b, am, bm, gamma), tspan, z0);
        z0 = z(3,:);
    end
end

% this function immitates the system
function dzdt = syst(time, z, ucom, a, b, am, bm, gamma)

    % this part of the code will change the parameters of they plant
    % during the simulation
    % change the condition from false to true to enable it
    if false
        if time > 200
            a = 3;
            b = 2;
        end
        
        if time > 400
            a = 4;
            b = 4;
        end
    end
    
    ym = z(1);
    yp = z(2);
    s = z(3);
    t0= z(4);
    
    u = t0* ucom - s * yp;
    e = yp - ym;
    
    dymdt= -am * ym + bm * ucom;
    dypdt= -a * yp + b * u;
    dsdt = gamma * e * b * yp;
    dtdt = -gamma * e * b * ucom;
    
    dzdt = [dymdt; dypdt; dsdt; dtdt];
end
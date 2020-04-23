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
gamma = 0.2;
ym0 = 0;
yp0 = 0;
z0 = [ym0; yp0; s0; t0];
%create time span and input command
t_sampling = 0.1;
times = 0:t_sampling:1200;

ucom = randn(length(times), 1);

[time, z] = ode45(@(time,z) syst(time, z, ucom, times, a, b, am, bm, gamma), times, z0);

% this function immitates the model
function dzdt = syst(time, z, ucomall, ucomt, a, b, am, bm, gamma)

    ucom = interp1(ucomt, ucomall, time);
    
    ym = z(1);
    yp = z(2);
    s = z(3);
    t = z(4);
    
    u = t * ucom - s * yp;
    e = yp - ym;
    
    dymdt = -am * ym + bm * ucom;
    dypdt = -a * yp + b * u;
    dsdt = gamma * e * b * yp;
    dtdt = -gamma * e * b * ucom;
    
    dzdt = [dymdt; dypdt; dsdt; dtdt];
end

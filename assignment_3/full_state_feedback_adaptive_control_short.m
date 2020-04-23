clear;

Ap = [0, 1; -16, -1];
Bp = [0; 1];

Am = [0, 1; -25, -10];
Q = 1;
Bm = Q*Bp;

syms p11 p12 p21 p22;
p = [p11, p12; p21, p22];
q = solve(Am'* p + p*Am == -3000*eye(2), {p11, p12, p21, p22});
p = [double(q.p11), double(q.p12); double(q.p21), double(q.p22)];

syms k1 k2;
K_star = [k1, k2];
q = solve(Ap + Bp*K_star == Am, {k1, k2});
K_star= [double(q.k1), double(q.k2)];


K0 = [0, 0];
xm0 = [0;0];
xp0 = [0;0];
z0 = [xm0; xp0; K0'];

%create time span and input command
t_sampling = 2;
time = 0:t_sampling:2000;

ucom = randn(length(time), 1);

[time, zs] = ode45(@(t, z) syst(t, z, ucom, time, Am, Bm, Ap, Bp, Q, p), time, z0);

xm1 = zs(:, 1);
xm2 = zs(:, 2);
xp1 = zs(:, 3);
xp2 = zs(:, 4);
K1 = zs(:, 5);
K2 = zs(:, 6);

figure;
subplot(3,1,1)
plot(time, [xp1-xm1]);
legend(["xm1", "xp1", "e"]);

subplot(3,1,2)
plot(time, [xp2-xm2]);
legend(["xm2", "xp2", "e"]);

subplot(3,1,3)
plot(time, [K1, K2]);
legend(["K1", "K2"]);


function dzdt = syst(t, z, ucomall, ucomt, Am, Bm, Ap, Bp, Q, p)
    xm = z(1:2);
    xp = z(3:4);
    K = z(5:6)';
    
    ucom = interp1(ucomt, ucomall, t);
    u = K * xp + Q * ucom;
    e = xp - xm;
    
    dxmdt = Am * xm + Bm * ucom;
    dxpdt = Ap * xp + Bp * u;
    dKdt = -Bp'*p*e*xp';
    
    dzdt = [dxmdt; dxpdt; dKdt'];
end
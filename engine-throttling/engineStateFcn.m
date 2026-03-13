

function xdot = engineStateFcn(x,u)
% x(1) = Ptc,  u(1) = mdot_O 
Ptc  = x(1);
mdot = u(1);

R = 282;
Tc = 3000; 
ao = 2.93e-4;
n = 0.518;
Vc = 9.5e-4;
r = 0.01;
Lf = 0.315;
rhof = 900;
gamma = 1.13;
nc = 0.7;
Cd = 1;
Ath = (pi*0.02^2)/4;
mdot = max(mdot,0);

K1 = (R*Tc)/Vc;
K2 = (2*pi*r*Lf*ao);
K3 = rhof;
K4 = R*Tc;

choked = (2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
K5 = (gamma*Cd*Ath)/(nc*sqrt(gamma*R*Tc)) * choked;

xdot = K1 * ( (K2/(pi*r^2)^n)*(mdot^n)*(K3 - (Ptc/K4)) + mdot - K5*Ptc );


end



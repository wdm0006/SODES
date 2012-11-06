% Hw 1 (2DOF)
% TR link
% Dynamics 
% Newton-Euler equations of motion

clear all
close all
clc

m1 = 1;
m2 = 1;
L  = 1;
IO = 1;
g = 9.81;

syms t p f21 

r = sym('r(t)');
theta = sym('theta(t)');

omega_ = [0 0 diff(theta,t)];
alpha_ = diff(omega_,t);

% Link 1
c =  cos(theta);
s =  sin(theta);
xP = p*c;
yP = p*s;
rP_ = [xP yP 0];
F21_ = [-f21*s f21*c 0];
eom1m_ = IO*alpha_ - cross(rP_,F21_);
eom1mz = eom1m_(3);

%Link 2 
xC = r*c;
yC = r*s;
rC_ = [xC yC 0];
vC_ = diff(rC_,t);
aC_ = diff(vC_,t);
G2_ = [0 -m2*g 0];
IC = m2*L^2/12;

eom2f_ = m2*aC_ + F21_ - G2_;
eom2fx = eom2f_(1);
eom2fy = eom2f_(2);
eom2m_ = IC*alpha_ - cross(rP_-rC_,-F21_);
eom2mz = eom2m_(3);

sol = solve(eom2fx, eom2mz, f21, p);
f21s = sol.f21;
ps = sol.p;

EqI  = subs(eom1mz,{'f21','p'},{f21s,ps});
EqII = subs(eom2fy,{'f21','p'},{f21s,ps});

% list for the symbolical variables 
slist={diff('theta(t)',t,2),diff('r(t)',t,2),...
       diff('theta(t)',t),diff('r',t),'theta(t)','r(t)'};
nlist = {'ddq1', 'ddq2', 'x(2)', 'x(4)', 'x(1)','x(3)'};

eq1 = subs(EqI,slist,nlist);
eq2 = subs(EqII,slist,nlist);

sole = solve(eq1,eq2,'ddq1, ddq2');
ddq1s = sole.ddq1;
ddq2s = sole.ddq2;
options= odeset('RelTol',1e-3,'MaxStep',1e-3);
statenames={'omega 1', 'theta 1', 'v', 'r'};
stateunits={'radians/sec', 'radians', 'meters/sec', 'meter'};
sodes({ddq1s,ddq2s},'x',0,5,[0.001,0,0.001,0], 'Integrator', 'ode45', 'Terminal', 'ON', 'Plotting', 'ON', 'StateNames', statenames, 'StateUnits', stateunits,'Options',options);

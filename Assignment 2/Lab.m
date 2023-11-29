clear all; close all; clc

Kq = 10;
a = 0.238;
g = 981;
K12 = 6;
Kh = 2;
A = 63.6;
h10 = 14;
h20 = 14;
Vsubs = h10/Kh;


syms u1 u2 V1 V2

f1 = Kq*u1/(A*Kh) - a*sqrt(2*g*Kh*V1)/(A*Kh) - K12*(Kh*V1-Kh*V2)/(A*Kh);
f2 = Kq*u2/(A*Kh) - a*sqrt(2*g*Kh*V2)/(A*Kh) + K12*(Kh*V1-Kh*V2)/(A*Kh);



f = [f1;f2];
 
V = [V1;V2];
u = [u1 ; u2];

deltaA = jacobian(f,V);
deltaA_V = simplify(subs(deltaA, V, [Vsubs;Vsubs]));
deltaB = double(jacobian(f,u));

%with respect to h
syms h1 h2

f1_h = Kq*u1/A - a*sqrt(2*g*h1)/A - K12*(h1-h2)/A;
f2_h = Kq*u2/A - a*sqrt(2*g*h2)/A + K12*(h1-h2)/A;
f_h = [f1_h; f2_h];
h = [h1;h2];
deltaA = jacobian(f_h,h);
deltaA_h = double(simplify(subs(deltaA, h, [h10;h10])));
deltaB = double(jacobian(f_h,u));
deltaC = eye(2);

f1_h_sim = - a*sqrt(2*g*h1)/A - K12*(h1-h2)/A;
f2_h_sim = - a*sqrt(2*g*h2)/A + K12*(h1-h2)/A;

%Stationary points 
u10 = solve(f1_h==0, u1);
u20 = solve(f2_h==0, u2);

u10 = double(subs(u10, h, [h10; h20]));
u20 = double(subs(u20, h, [h10; h20]));
%% controllability and 

S =  ctrb(deltaA_h, deltaB);
Obs = obsv(deltaA_h, deltaC);

rank(S)
rank(Obs)

model = 'Simulink_lab1';
model_data = sim(model);
sim_time = model_data.simout.time;
sim_data = model_data.simout.signals.values;
t_end= find(sim_data(:,1) > 14.137592,1,'first');
tau = (t_end - 101)*0.1;

%% SISO-system

f1_h = Kq*u1/A - a*sqrt(2*g*h1)/A;
f2_h = Kq*u2/A - a*sqrt(2*g*h2)/A;
f_h = [f1_h; f2_h];
h = [h1;h2];
deltaA = jacobian(f_h,h);
deltaA_SISO= double(simplify(subs(deltaA, h, [h10;h10])));
deltaB_SISO = double(jacobian(f_h,u));

statespace = ss(deltaA_SISO, deltaB_SISO, deltaC,0);

trans = tf(statespace);
%% PI-controller

lambda = tau/10;
tau2 = trans(1, 1).Denominator{1, 1}(1)/trans(1, 1).Denominator{1, 1}(2);
Kp = trans(1, 1).Numerator{1}(2)/ trans(1, 1).Denominator{1, 1}(2);
Ki =  0.8;   %Kp/tau2;


%% LQ controller 

rho = 10;
Qu =  eye(size(deltaA_h,1));
Qx = rho*Qu;
%L = inv(Qu)*(deltaB'*S);
[Lc,K,P] = lqr(ss(deltaA_h, deltaB,deltaC,0), Qx, Qu); 
Kr = -inv(inv(deltaA_h-deltaB*Lc)*deltaB);

%% LQI-controller
eta = 0.1;
Atot = [deltaA_h zeros(2);-deltaC zeros(2)];
Btot = [deltaB;zeros(2)];
Qtot = blkdiag(rho*eye(2),eta*eye(2));
Ltot= lqr(Atot, Btot, Qtot, Qu);
Lc = Ltot(:,1:2)
Li = Ltot(:,3:4)

Kr  = [Lc eye(2)]*inv([deltaA_h deltaB;deltaC eye(2)])*[zeros(3,1);1]

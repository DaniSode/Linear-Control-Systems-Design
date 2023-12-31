clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              %
%                LAB Assignment                %
%    Linear Control Systems Design - SSY285    %
%                   Authors:                   %
%              Daniel Söderqvist               %
%               Johannes Lundahl               %
%             Martienn Sigthorsson             %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Exercise 3

% Definitions of symbols and constants
syms h_1 h_2 u_1 u_2 V_1 V_2 real
At = 63.6;
K_q = 10;
a = 0.238;
g = 981;
K_12 = 6;
K_h = 2;

% Define state equations
disp('Numerical equations for simulink:')
eqs = [(1/(At))*(K_q*u_1 - a*sqrt(2*g*h_1) - K_12*(h_1 - h_2));
       (1/At)*(K_q*u_2 - a*sqrt(2*g*h_2) + K_12*(h_1 - h_2))]

% Define state vector and input vector
x = [h_1; h_2];
u = [u_1; u_2];

% Linearize the system by finding Jacobian matrix
A_lin = jacobian(eqs, x);
B_lin = jacobian(eqs, u);


% Linearized state equations
disp('A matrix:'); disp(A_lin)
disp('B matrix:'); disp(B_lin)

% Given linearization point
point = [14, 14];

% Substitute
A_lin_sub = subs(A_lin, x', point);
B_lin_sub = subs(B_lin, x', point);

% Linearized A and B matrices with substitution
disp('Linearized A matrix'); disp(double(A_lin_sub))
disp('Linearized B matrix'); disp(double(B_lin_sub))

% Solve for init
eqs_solve = subs(eqs, x', point);
u_1_init = double(solve(eqs_solve(1)==0,u_1));
u_2_init = double(solve(eqs_solve(2)==0,u_2));

% Extras for lambda tuning
eqs_simple = [(1/(At))*(K_q*u_1 - a*sqrt(2*g*h_1));
           (1/At)*(K_q*u_2 - a*sqrt(2*g*h_2))];
A_simple = double(subs(jacobian(eqs_simple, x), x', [14, 14]));
B_simple = double(subs(jacobian(eqs_simple, u), x', [14, 14]));

%% Exercise 4

% Define state equations
eqs_v = [(1/(At*K_h))*(K_q*u_1 - a*sqrt(2*g*K_h*V_1) - K_12*(K_h*V_1 - K_h*V_2));
       (1/(At*K_h))*(K_q*u_2 - a*sqrt(2*g*K_h*V_2) + K_12*(K_h*V_1 - K_h*V_2))];

% Define state vector
x_v = [V_1; V_2];

% Linearize the system by finding Jacobian matrix
A_lin_v = jacobian(eqs_v, x_v);
B_lin_v = jacobian(eqs_v, u);

% Linearized state equations
disp('A matrix with voltage:'); disp(A_lin_v)
disp('B matrix with voltage:'); disp(B_lin_v)
disp('T matrix:'); disp(B_lin*B_lin_v^(-1))

%% Exercise 5

% Matrices for state space model
disp('Matrices:')
A = double(A_lin_sub)
B = double(B_lin_sub)
C = eye(2)
D = zeros(2)

obs = obsv(A,C);
ctr = ctrb(A,B);
disp('Observability matrix:'); disp(obs)
disp('Controllability matrix:'); disp(ctr)
disp('Rank of Observability matrix:'); disp(rank(obs))
disp('Rank of Controllability matrix:'); disp(rank(ctr))

%% Exercise 6
ts = 0.1;
sysc = ss(A,B,C,D);
sysd = c2d(sysc, ts);

Simulation_Time = 300;
open("controlstructureforassignment.mdl")
sim_model = sim("Simulation.slx", 'ReturnWorkspaceOutputs', 'on');

%% Exercise 7
disp('Transfer function for the system:')
sys_simple = ss(A_simple,B_simple,C,D)
transferfcn = tf(sys_simple)

h10 = 14;
h20 = 14;
u10 = u_1_init;
u20 = u_2_init;

% Decide these parameters
% PI
Tau = 1/0.02215;
lambda = Tau/10;
Kp = 0.1572*Tau;
Ki = Tau/(Kp*lambda);
Ki = Kp/Tau;
% Finetuning
%Kp = 10;
%Ki = 1.1;


%% Exercise 8 and 9

rho = 10;
Qu =  eye(size(A,1));
Qx = rho*Qu;
[Lc,K,P] = lqr(ss(A,B,C,D), Qx, Qu); 
Kr = -pinv(pinv(A-B*Lc)*B);

%% Exercise 10 

eta = 0.1;
Atot = [A,zeros(2);-C,zeros(2)];
Btot = [B; zeros(2)];
Qtot = blkdiag(rho*eye(2),eta*eye(2));
[Lctot,Ktot,Ptot] = lqr(Atot,Btot,Qtot,Qu);
L = Lctot(:,1:2);
Li = Lctot(:,3:4);
Kur = -pinv(C*pinv(A)*B);
Kxr = pinv(A)*B*pinv((C*pinv(A)*B));



clear all; close all; clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              %
%              Home Assignment 3               %
%    Linear Control Systems Design - SSY285    %
%                   Authors:                   %
%              Daniel Söderqvist               %
%               Johannes Lundahl               %
%              Martienn Sigborsson             %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Definitions

% Definitions of variables as syms
syms i_a v_a R L K_E K_T omega_1 omega_2 omega_3 phi_1 phi_2 phi_3 J_1 J_2 D_1 D_2 B_f T_e real positive
syms i_a_dot omega_1_dot omega_2_dot omega_3_dot phi_1_dot phi_2_dot phi_3_dot real positive
T_a = K_T*i_a;
e_a = K_E*omega_1;

% Define state and input vector
x = [i_a; phi_1; phi_2; phi_3; omega_1; omega_2];
xp = [i_a_dot; phi_1_dot; phi_2_dot; phi_3_dot; omega_1_dot; omega_2_dot];
u = [v_a; T_e];

% Define system equations
equations = [v_a == R*i_a + L*i_a_dot + e_a;
             J_1*omega_1_dot == T_a - D_1*(phi_1 - phi_2);
             J_2*omega_2_dot == D_1*(phi_1 - phi_2) - D_2*(phi_2 - phi_3);
             phi_1_dot == omega_1;
             phi_2_dot == omega_2;
             phi_3_dot == (D_2*(phi_2 - phi_3) + T_e)/B_f]; 

equations = lhs(equations) - rhs(equations);

% Define A and B matrices for hole system
Am = -jacobian(equations, xp) \ jacobian(equations, x);
Bm = -jacobian(equations, xp) \ jacobian(equations, u);

% Neglect term with L solve for i_a and insert for all other i_a's
first_eq = subs(equations(1), L, 0);
i_a_new = solve(first_eq, i_a);
equations_new = subs(equations(2:end), i_a, i_a_new);

% Define new A and B matrices neglecting L term
A = -jacobian(equations_new, xp(2:end)) \ jacobian(equations_new, x(2:end));
B = -jacobian(equations_new, xp(2:end)) \ jacobian(equations_new, u);
% Define C and D matrices for first case
y_1 = [phi_2; omega_2];
C_a = jacobian(x(2:end), y_1)';
D_a = jacobian(u, y_1);

% Define C and D matrices for second case
y_2 = [i_a; phi_3_dot];
C_b = [-jacobian(first_eq, y_2(1)) \ jacobian(first_eq, x(2:end));
       -jacobian(equations_new(5), y_2(2)) \ jacobian(equations_new(5), x(2:end))];
D_b = [-jacobian(first_eq, y_2(1)) \ jacobian(first_eq, u);
       -jacobian(equations_new(5), y_2(2)) \ jacobian(equations_new(5), u)];


%% Substitution
variables = [K_E; K_T; J_1; J_2; B_f; D_2];
values = [10^(-1); 10^(-1); 10^(-5); 4*10^(-5); 2*10^(-3); 2];

% Display new matrices with numerical values
A_new = subs(A, variables, values);
B_new = subs(B, variables, values);
C_new_a = subs(C_a, variables, values);
D_new_a = subs(D_a, variables, values);
C_new_b = subs(C_b, variables, values);
D_new_b = subs(D_b, variables, values);

variables_c = [R; D_1];
values_c = [1; 20];

% Display new matrices with numerical values
A_c = double(subs(A_new, variables_c, values_c));
B_c = double(subs(B_new, variables_c, values_c));
C_c_1 = double(subs(C_new_a, variables_c, values_c));
D_c_1 = double(subs(D_new_a, variables_c, values_c));
C_c_2 = double(subs(C_new_b, variables_c, values_c));
D_c_2 = double(subs(D_new_b, variables_c, values_c));


% For controlling 
h = 0.01;
sys = ss(A_c, B_c, C_c_1, D_c_1);
sys_d = c2d(sys, h);
[Ad,Bd,C,D] = ssdata(sys_d)

% Augmented system 
C_omega = C(2,:);
Ad_aug = [Ad zeros(5,1); -C_omega 1];
Bd_aug = [Bd; zeros(1,2)];
Cd_aug = [C_omega 0];



% Time discrete white noise
Qw = [0.01 0 ; 0 1.111e-3];
Qv = [4.44e-5 0; 0 1.11e-5];
Nn = 0; 


new_sysd = ss(Ad,[Bd Bd],C,0,h)

% Kalman filter
[~, K, P] = kalman(new_sysd, Qw, Qv,Nn);


observer_eigin = eig(Ad-K*C)




Qx = diag([ 1, 10, 1, 1/1000, 1/1000, 11]);
Qu = 100000*diag([1, 1]);








% LQ solution
Ltot = dlqr(Ad_aug, Bd_aug,Qx, Qu);

L = Ltot(:,1:5)
Li = Ltot(:,6)









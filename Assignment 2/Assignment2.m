clear all; close all; clc

%% Definitions

% Definitions of variables as syms
syms i_a v_a_2 R L K_E K_T omega_1 omega_2 omega_3 phi_1 phi_2 phi_3 J_1 J_2 D_1 D_2 B_f T_e_2 real
syms i_a_dot omega_1_dot omega_2_dot omega_3_dot phi_1_dot phi_2_dot phi_3_dot real
T_a = K_T*i_a;
e_a = K_E*omega_1;

% Define state and input vector
x = [i_a; phi_1; omega_1; phi_2; omega_2; phi_3];
xp = [i_a_dot; phi_1_dot; omega_1_dot; phi_2_dot; omega_2_dot; phi_3_dot];
u_2 = [v_a_2; T_e_2];

% Define system equations
equations = [v_a_2 == R*i_a + L*i_a_dot + e_a;
             J_1*omega_1_dot == T_a - D_1*(phi_1 - phi_2);
             J_2*omega_2_dot == D_1*(phi_1 - phi_2) - D_2*(phi_2 - phi_3);
             phi_1_dot == omega_1;
             phi_2_dot == omega_2;
             phi_3_dot == (D_2*(phi_2 - phi_3) + T_e_2)/B_f]; 
disp('Equations:')
equations = lhs(equations) - rhs(equations)

% Define A and B matrices for hole system
disp('All considered:')
Am = -jacobian(equations, xp) \ jacobian(equations, x)
Bm = -jacobian(equations, xp) \ jacobian(equations, u_2)

% Neglect term with L solve for i_a and insert for all other i_a's
first_eq = subs(equations(1), L, 0);
i_a_new = solve(first_eq, i_a);
equations_new = subs(equations(2:end), i_a, i_a_new);

% Define new A and B matrices neglecting L term
disp('Considering L = 0:')
A = -jacobian(equations_new, xp(2:end)) \ jacobian(equations_new, x(2:end))
B = -jacobian(equations_new, xp(2:end)) \ jacobian(equations_new, u_2)

% Define C and D matrices for first case
disp('Case 1:')
y_1 = [phi_2; omega_2];
C_a = jacobian(x(2:end), y_1)'
D_a = jacobian(u_2, y_1)

% Define C and D matrices for second case
disp('Case 2:')
y_2 = [i_a; phi_3_dot];
C_b = [-jacobian(first_eq, y_2(1)) \ jacobian(first_eq, x(2:end));
       -jacobian(equations_new(5), y_2(2)) \ jacobian(equations_new(5), x(2:end))]
D_b = [-jacobian(first_eq, y_2(1)) \ jacobian(first_eq, u_2);
       -jacobian(equations_new(5), y_2(2)) \ jacobian(equations_new(5), u_2)]

% Simulation case
C_sim = [jacobian(x(2:end), [omega_1, omega_2])';
         -jacobian(equations_new(5), phi_3_dot) \ jacobian(equations_new(5), x(2:end))];
D_sim = [jacobian(u_2, y_1);
         -jacobian(equations_new(5), phi_3_dot) \ jacobian(equations_new(5), u_2)];

%% Substitution
variables = [K_E; K_T; J_1; J_2; B_f; D_2];
values = [10^(-1); 10^(-1); 10^(-5); 4*10^(-5); 2*10^(-3); 2];

% Display new matrices with numerical values
disp('A with numerical values:')
A_new = subs(A, variables, values);
B_new = subs(B, variables, values);
C_new_a = subs(C_a, variables, values);
D_new_a = subs(D_a, variables, values);
C_new_b = subs(C_b, variables, values);
D_new_b = subs(D_b, variables, values);

%% Controllability 
syms lam real
Ilambda = eye(size(A,1))*lam;
Contr = [A-Ilambda, B];

%Case 1
Ctr_1 = [B_new, A_new*B_new, A_new^2*B_new, A_new^3*B_new, A_new^4*B_new];

ranken = rank(Ctr_1);
echolon = rref(Ctr_1);



% if rank(Ctr_1)~= 0
%     fprintf('System with controllability matrix Ctr_1 has rank %g and is controllable!', rank(Ctr_1))
% else 
%     fprintf('System with controllability matrix Ctr_1 has rank %g and is NOT controllable', rank(Ctr_1))
% end


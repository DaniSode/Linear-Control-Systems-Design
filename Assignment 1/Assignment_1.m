clear all; close all; clc

syms i_a v_a R L K_E K_T omega_1 omega_2 omega_3 phi_1 phi_2 phi_3 J_1 J_2 D_1 D_2 B_f T_e real
syms i_a_dot omega_1_dot omega_2_dot omega_3_dot phi_1_dot phi_2_dot phi_3_dot real
T_a = K_T*i_a;
e_a = K_E*omega_1;

x = [i_a; phi_1; omega_1; phi_2; omega_2; phi_3];
xp = [i_a_dot; phi_1_dot; omega_1_dot; phi_2_dot; omega_2_dot; phi_3_dot];
u = [v_a; T_e];

equations = [v_a == R*i_a + L*i_a_dot + e_a;
             J_1*omega_1_dot == T_a - D_1*(phi_1 - phi_2);
             J_2*omega_2_dot == D_1*(phi_1 - phi_2) - D_2*(phi_2 - phi_3);
             phi_1_dot == omega_1;
             phi_2_dot == omega_2;
             phi_3_dot == (D_2*(phi_2 - phi_3) + T_e)/B_f]; 

disp('Equations:')
equations = lhs(equations) - rhs(equations)

disp('All considered:')
Am = -jacobian(equations, xp) \ jacobian(equations, x)
Bm = -jacobian(equations, xp) \ jacobian(equations, u)

% Without L term
first_eq = subs(equations(1), L, 0);
i_a_new = solve(first_eq, i_a);
equations_new = subs(equations(2:end), i_a, i_a_new);

disp('Considering L = 0:')
A = -jacobian(equations_new, xp(2:end)) \ jacobian(equations_new, x(2:end))
B = -jacobian(equations_new, xp(2:end)) \ jacobian(equations_new, u)

disp('Case 1:')
y_1 = [phi_2; omega_2];
C_a = jacobian(x(2:end), y_1)'
D_a = jacobian(u, y_1)

% We solve for i_a and we know that omega_3 = phi_3_dot
disp('Case 2:')
y_2 = [i_a; phi_3_dot];
C_b = [-jacobian(first_eq, y_2(1)) \ jacobian(first_eq, x(2:end));
       -jacobian(equations_new(5), y_2(2)) \ jacobian(equations_new(5), x(2:end))]
D_b = [-jacobian(first_eq, y_2(1)) \ jacobian(first_eq, u);
       -jacobian(equations_new(5), y_2(2)) \ jacobian(equations_new(5), u)]

variables = [R; K_E; K_T; J_1; J_2; B_f; D_1; D_2];
values = [1; 10^(-1); 10^(-1); 10^(-5); 4*10^(-5); 2*10^(-3); 20; 2];

disp('A with numerical values:')
A_new = double(subs(A, variables, values))
B_new = double(subs(B, variables, values));
C_new_a = double(subs(C_a, variables, values));
D_new_a = double(subs(D_a, variables, values));
C_new_b = double(subs(C_b, variables, values));
D_new_b = double(subs(D_b, variables, values));

disp('Eigenvalues:')
double(eig(A_new))

ss_model_1 = ss(A_new, B_new, C_new_a, D_new_a);
ss_model_2 = ss(A_new, B_new, C_new_b, D_new_b);

subplot(2,1,1)
pzplot(ss_model_1)

subplot(2,1,2)
pzplot(ss_model_2)
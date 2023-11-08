clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              %
%              Home Assignment 1               %
%    Linear Control Systems Design - SSY285    %
%                   Authors:                   %
%              Daniel Söderqvist               %
%               Johannes Lundahl               %
%              Martienn Sigborsson             %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Definitions

% Definitions of variables as syms
syms i_a v_a R L K_E K_T omega_1 omega_2 omega_3 phi_1 phi_2 phi_3 J_1 J_2 D_1 D_2 B_f T_e real
syms i_a_dot omega_1_dot omega_2_dot omega_3_dot phi_1_dot phi_2_dot phi_3_dot real
T_a = K_T*i_a;
e_a = K_E*omega_1;

% Define state and input vector
x = [i_a; phi_1; omega_1; phi_2; omega_2; phi_3];
xp = [i_a_dot; phi_1_dot; omega_1_dot; phi_2_dot; omega_2_dot; phi_3_dot];
u = [v_a; T_e];

% Define system equations
equations = [v_a == R*i_a + L*i_a_dot + e_a;
             J_1*omega_1_dot == T_a - D_1*(phi_1 - phi_2);
             J_2*omega_2_dot == D_1*(phi_1 - phi_2) - D_2*(phi_2 - phi_3);
             phi_1_dot == omega_1;
             phi_2_dot == omega_2;
             phi_3_dot == (D_2*(phi_2 - phi_3) + T_e)/B_f]; 
disp('Equations:')
equations = lhs(equations) - rhs(equations)

% Define A and B matrices for hole system
disp('All considered:')
Am = -jacobian(equations, xp) \ jacobian(equations, x)
Bm = -jacobian(equations, xp) \ jacobian(equations, u)

% Neglect term with L solve for i_a and insert for all other i_a's
first_eq = subs(equations(1), L, 0);
i_a_new = solve(first_eq, i_a);
equations_new = subs(equations(2:end), i_a, i_a_new);

% Define new A and B matrices neglecting L term
disp('Considering L = 0:')
A = -jacobian(equations_new, xp(2:end)) \ jacobian(equations_new, x(2:end))
B = -jacobian(equations_new, xp(2:end)) \ jacobian(equations_new, u)

% Define C and D matrices for first case
disp('Case 1:')
y_1 = [phi_2; omega_2];
C_a = jacobian(x(2:end), y_1)'
D_a = jacobian(u, y_1)

% Define C and D matrices for second case
disp('Case 2:')
y_2 = [i_a; phi_3_dot];
C_b = [-jacobian(first_eq, y_2(1)) \ jacobian(first_eq, x(2:end));
       -jacobian(equations_new(5), y_2(2)) \ jacobian(equations_new(5), x(2:end))]
D_b = [-jacobian(first_eq, y_2(1)) \ jacobian(first_eq, u);
       -jacobian(equations_new(5), y_2(2)) \ jacobian(equations_new(5), u)]


%% Substitution
variables = [R; K_E; K_T; J_1; J_2; B_f; D_1; D_2];
values = [1; 10^(-1); 10^(-1); 10^(-5); 4*10^(-5); 2*10^(-3); 20; 2];

% Display new matrices with numerical values
disp('A with numerical values:')
A_new = double(subs(A, variables, values))
B_new = double(subs(B, variables, values));
C_new_a = double(subs(C_a, variables, values));
D_new_a = double(subs(D_a, variables, values));
C_new_b = double(subs(C_b, variables, values));
D_new_b = double(subs(D_b, variables, values));

% Calculate and show eigenvalues
disp('Eigenvalues:')
eigenvalues = double(eig(A_new))

% Define state-space models of the 2 cases
ss_model_1 = ss(A_new, B_new, C_new_a, D_new_a);
ss_model_2 = ss(A_new, B_new, C_new_b, D_new_b);

% Plot the poles for the cases
subplot(2,1,1)
pzplot(ss_model_1)
subplot(2,1,2)
pzplot(ss_model_2)


%% Define transfer function
s = tf("s");
G_s = C_new_b*inv(s*eye(size(A_new,1))-A_new)*B_new + D_new_b

% Display poles and transmission zeros
format longg
disp('Poles:')
poles = pole(ss_model_2)
disp('Transmission zeros:')
trans_zeros = tzero(ss_model_2)


%% Simulatation
time = linspace(0, 0.05, 101);
v_a = [0*ones(1, (length(time)-1)*0.05), 10*ones(1, length(time)-(length(time)-1)*0.05)];
T_e = [0*ones(1, (length(time)-1)*0.5), -0.1*ones(1, length(time)-(length(time)-1)*0.5)];
u = [v_a; T_e];

% Simulation of the 2 models
[y_1, t_1, x_1] = lsim(ss_model_1, u, time);
[y_2, t_2, x_2] = lsim(ss_model_2, u, time);

% Plotting
figure
subplot(3,2,1)
plot(t_2, x_1(:, 2))
axis padded
title('Omega 1'); xlabel('Time [s]'); ylabel('rad/s')
subplot(3,2,3)
plot(t_2, x_1(:, 4))
axis padded
title('Omega 2'); xlabel('Time [s]'); ylabel('rad/s')
subplot(3,2,5)
plot(t_1, y_2(:, 2))
axis padded
title('Omega 3'); xlabel('Time [s]'); ylabel('rad/s')
subplot(3,2,2)
plot(t_1, y_1(:, 1))
axis padded
title('Current'); xlabel('Time [s]'); ylabel('i')
subplot(3,2,4)
plot(time, u(1, :))
axis padded
title('Input voltage'); xlabel('Time [s]'); ylabel('v')
subplot(3,2,6)
plot(time, u(2, :))
axis padded
title('External forces'); xlabel('Time [s]'); ylabel('Nm')



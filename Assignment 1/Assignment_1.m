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
C_new_sim = double(subs(C_sim, variables, values));
D_new_sim = double(subs(D_sim, variables, values));


% Calculate and show eigenvalues
disp('Eigenvalues of A:')
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
disp('Transfer function for second model:')
sys_1 = tf(ss_model_1);
sys_2 = tf(ss_model_2)

% Display poles and transmission zeros
format longg
disp('Poles for second model:')
poles = pole(ss_model_2)
disp('Transmission zeros for second model:')
trans_zeros = tzero(ss_model_2)


%% Simulatation
% Creating model and testing to simulate when ramping input
ss_model_sim = ss(A_new, B_new, C_new_sim, D_new_sim);
t_1 = 0:0.001:0.5;
T_e_1 = min(0,max(-7*(t_1-0.25),-0.1));
v_a_1 = max(0,min(600*(t_1-0.025),10));
u_1 = [v_a_1; T_e_1];
[y_2, t_sim_2] = lsim(ss_model_sim, u_1, t_1);

% Simulating when stepwise changing input
t_2 = linspace(0, 0.05, 1001);
v_a_2 = [0*ones(1, (length(t_2)-1)*0.05), 10*ones(1, length(t_2)-(length(t_2)-1)*0.05)];
T_e_2 = [0*ones(1, (length(t_2)-1)*0.5), -0.1*ones(1, length(t_2)-(length(t_2)-1)*0.5)];
u_2 = [v_a_2; T_e_2];
[y_1, t_sim_1] = lsim(ss_model_sim, u_2, t_2);

% Plotting
figure
subplot(4,3,1)
plot(t_sim_2, y_2(:, 1))
axis padded; grid on
title('Omega 1 (ramped input)'); xlabel('Time [s]'); ylabel('rad/s')
subplot(4,3,2)
plot(t_sim_2, y_2(:, 2))
axis padded; grid on
title('Omega 2 (ramped input)'); xlabel('Time [s]'); ylabel('rad/s')
subplot(4,3,3)
plot(t_sim_2, y_2(:, 3))
axis padded; grid on
title('Omega 3 (ramped input)'); xlabel('Time [s]'); ylabel('rad/s')
subplot(4,3,7)
plot(t_sim_1, y_1(:, 1))
axis padded; grid on
title('Omega 1 (step input)'); xlabel('Time [s]'); ylabel('rad/s')
subplot(4,3,8)
plot(t_sim_1, y_1(:, 2))
axis padded; grid on
title('Omega 2 (step input)'); xlabel('Time [s]'); ylabel('rad/s')
subplot(4,3,9)
plot(t_sim_1, y_1(:, 3))
axis padded; grid on
title('Omega 3 (step input)'); xlabel('Time [s]'); ylabel('rad/s')
subplot(4,2,3)
plot(t_1, u_1(1, :))
axis padded; grid on
title('Input voltage'); xlabel('Time [s]'); ylabel('v')
subplot(4,2,4)
plot(t_1, u_1(2, :))
axis padded; grid on
title('External forces'); xlabel('Time [s]'); ylabel('Nm')
subplot(4,2,7)
plot(t_2, u_2(1, :))
axis padded; grid on
title('Input voltage'); xlabel('Time [s]'); ylabel('v')
subplot(4,2,8)
plot(t_2, u_2(2, :))
axis padded; grid on
title('External forces'); xlabel('Time [s]'); ylabel('Nm')



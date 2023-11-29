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

% Simulation case
C_sim = [-jacobian(first_eq, i_a) \ jacobian(first_eq, x(2:end));
         jacobian(x(2:end), [phi_1, phi_2, phi_3])';
         jacobian(x(2:end), [omega_1, omega_2])';
         -jacobian(equations_new(5), phi_3_dot) \ jacobian(equations_new(5), x(2:end))];
D_sim = [-jacobian(first_eq, i_a) \ jacobian(first_eq, u);
         jacobian(u, [phi_1, phi_2, phi_3])';
         jacobian(u, [omega_1, omega_2])';
         -jacobian(equations_new(5), phi_3_dot) \ jacobian(equations_new(5), u)];

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

% Define state-space models of the 2 cases without external forces
A_mod = A_new;
B_mod = -jacobian(equations_new, xp(2:end)) \ jacobian(equations_new, u(1));
B_mod = double(subs(B_mod, variables, values));
C_mod_a = C_new_a; 
C_mod_b = C_new_b;
D_mod_a = jacobian(u(1), y_1);
D_mod_a = double(subs(D_mod_a, variables, values))';
D_mod_b = [-jacobian(first_eq, y_2(1)) \ jacobian(first_eq, u(1));
           -jacobian(equations_new(5), y_2(2)) \ jacobian(equations_new(5), u(1))];
D_mod_b = double(subs(D_mod_b, variables, values));
ss_model_3 = ss(A_mod, B_mod, C_mod_a, D_mod_a);
ss_model_4 = ss(A_mod, B_mod, C_mod_b, D_mod_b);

% Plot the poles for the cases
figure
pzplot(ss_model_1(1,1),'g' ,ss_model_1(1,2),'b',ss_model_1(2,1),'k',ss_model_1(2,2),'b')
xlim([-1050 50])
ylim([-1650 1650])
title('Pole-Zero Map Case 1')
saveas(gcf, 'case1', 'epsc')
figure
pzplot(ss_model_2(1,1),'g' ,ss_model_2(1,2),'r',ss_model_2(2,1),'k',ss_model_2(2,2),'b')
xlim([-1000 50])
ylim([-1750 1750])
title('Pole-Zero Map Case 2')
saveas(gcf, 'case2', 'epsc')

% Plot the poles for the cases without external forces
figure
pzplot(ss_model_3(1,1),'b' ,ss_model_3(2,1),'r')
xlim([-1050 50])
ylim([-1650 1650])
title('Pole-Zero Map Case 1 without external forces')
saveas(gcf, 'case1ex', 'epsc')
figure
pzplot(ss_model_4(1,1),'b' ,ss_model_4(2,1),'r')
xlim([-1000 50])
ylim([-1750 1750])
title('Pole-Zero Map Case 2 without external forces')
saveas(gcf, 'case2ex', 'epsc')

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
% Creating model and simulatating
ss_model_sim = ss(A_new, B_new, C_new_sim, D_new_sim);
t = 0:0.0001:0.05;
T_e_sim = min(0,max(-0.1/0.0001*(t-0.025),-0.1));
v_a_sim = max(0,min(10/0.0001*(t-0.0025),10));
u_sim = [v_a_sim; T_e_sim];
[y, t_sim] = lsim(ss_model_sim, u_sim, t);

% Plotting
figure
subplot(3,3,1)
plot(t_sim, y(:, 2))
axis padded; grid on
title('Phi 1'); xlabel('Time [s]'); ylabel('rad')
subplot(3,3,2)
plot(t_sim, y(:, 3))
axis padded; grid on
title('Phi 2'); xlabel('Time [s]'); ylabel('rad')
subplot(3,3,3)
plot(t_sim, y(:, 4))
axis padded; grid on
title('Phi 3'); xlabel('Time [s]'); ylabel('rad')
subplot(3,3,4)
plot(t_sim, y(:, 5))
axis padded; grid on
title('Omega 1'); xlabel('Time [s]'); ylabel('rad/s')
subplot(3,3,5)
plot(t_sim, y(:, 6))
axis padded; grid on
title('Omega 2'); xlabel('Time [s]'); ylabel('rad/s')
subplot(3,3,6)
plot(t_sim, y(:, 7))
axis padded; grid on
title('Omega 3'); xlabel('Time [s]'); ylabel('rad/s')
subplot(3,3,7)
plot(t_sim, y(:, 1))
axis padded; grid on
title('Current'); xlabel('Time [s]'); ylabel('i')
subplot(3,3,8)
plot(t, u_sim(1, :))
axis padded; grid on
title('Input voltage'); xlabel('Time [s]'); ylabel('v')
subplot(3,3,9)
plot(t, u_sim(2, :))
axis padded; grid on
title('External forces'); xlabel('Time [s]'); ylabel('Nm')
saveas(gcf, 'sim_fig', 'epsc')

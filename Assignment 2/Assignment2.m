clear all; close all; clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              %
%              Home Assignment 2               %
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

equations = lhs(equations) - rhs(equations);

% Define A and B matrices for hole system
disp('System taken from MA1:')
Am = -jacobian(equations, xp) \ jacobian(equations, x);
Bm = -jacobian(equations, xp) \ jacobian(equations, u);

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
variables = [K_E; K_T; J_1; J_2; B_f; D_2];
values = [10^(-1); 10^(-1); 10^(-5); 4*10^(-5); 2*10^(-3); 2];

% Display new matrices with numerical values
A_new = subs(A, variables, values);
B_new = subs(B, variables, values);
C_new_a = subs(C_a, variables, values);
D_new_a = subs(D_a, variables, values);
C_new_b = subs(C_b, variables, values);
D_new_b = subs(D_b, variables, values);

%% a)
syms lam real positive
Ilambda = eye(size(A,1))*lam;
Contr = [A-Ilambda, B];

%CASE 1
%Controllability
S1 = [B_new, A_new*B_new, A_new^2*B_new, A_new^3*B_new, A_new^4*B_new];
%Observability 
O1 = [C_new_a; C_new_a*A_new; C_new_a*A_new^2; C_new_a*A_new^3; ...
    C_new_a*A_new^4];
disp('Case 1:')
disp('Rank of controllability:')
rankS1 = rank(S1)
disp('Rank of observability:')
rankO1 = rank(O1)
disp('Echelon matrix for controllability:')
R_S1 = rref(S1)
disp('Echelon matrix for observability:')
R_O1 = rref(O1)

%CASE 2
%Controllability
S2 = [B_new, A_new*B_new, A_new^2*B_new, A_new^3*B_new, A_new^4*B_new];
%Observability 
O2 = [C_new_b; C_new_b*A_new; C_new_b*A_new^2; C_new_b*A_new^3; ...
    C_new_b*A_new^4];
disp('Case 2:')
disp('Rank of controllability:')
rankS2 = rank(S2)
disp('Rank of observability:')
rankO2 = rank(O2)
disp('Echelon matrix for controllability:')
R_S2 = rref(S2)
disp('Echelon matrix for observability:')
R_O2 = rref(O2)


%% b)
% If the system is controllable, i.e. ctr1 has full rank it implies that 
% the system is also stabalizable
% Knowing fewer measurements ie lower observability rank means that we
% still can estimate states thats not observable if the system is still
% full rank when it comes to controllability. 

% Full rank in controllabiltiy means its fully stabilizable
% Full observabiltiy means fully detectable
% However in case 2 its not full rank which means we have to chekc if its
% detectable. 
disp('Eigenvalues of A:')
eigen = simplify(eig(simplify(A_new)))
svd([A_new;C_new_b])        % PBH test with lambda=0, results in rank<=4 (rank deficient)

%% C

variables_c = [R; D_1];
values_c = [1; 20];

% Display new matrices with numerical values
A_c = subs(A_new, variables_c, values_c);
B_c = subs(B_new, variables_c, values_c);
C_c_1 = subs(C_new_a, variables_c, values_c);
D_c_1 = subs(D_new_a, variables_c, values_c);
C_c_2 = subs(C_new_b, variables_c, values_c);
D_c_2 = subs(D_new_b, variables_c, values_c);

%CASE 1
%Controllability using built in functions 
disp('Rank of controllability after substituting both cases:')
S_c = ctrb(double(A_c),double(B_c));
rank_c1 = rank(S_c)
%Observability using built in functions 
disp('Rank of observability after substituting for case 1:')
O_c_1 = obsv(double(A_c), double(C_c_1));
rank_oc1 = rank(O_c_1)
%CASE 2
%Observability using built in functions
disp('Rank of observability after substituting for case 2:')
O_c_2 = obsv(double(A_c), double(C_c_2));
rank_oc2 = rank(O_c_2)

%% d) (Question: is it allowed to use c2d and compute the discrete time system matrix?)
syms s
Ts = 0.001;
A_d = expm(double(A_c)*Ts)

sys = ss(double(A_c), double(B_c), double(C_c_1), double(D_c_1));

sys_d = c2d(sys, Ts);
[sys_d_A,~,~,~] = ssdata(sys_d)

%d)
sysr = minreal(sys)


% ENPM 667 Fall 2022
% Final Project
% Jerry Pittman (117077120
% Bob Reiter (119457924)

clear; 
clc;
close all;

%------First Component------------
% Part A: Equations of motion
syms m1 m2 l1 l2 M g t F

syms x(t) theta1(t) theta2(t)

dx = diff(x);
dth1 = diff(theta1);
dth2 = diff(theta2);

%---Solve for Kinetic Energy
v1_sq = (dx + l1 * cos(theta1) * dth1).^2 + (l1 * sin(theta1) * dth1).^2;
v2_sq = (dx + l2 * cos(theta2) * dth2).^2 + (l2 * sin(theta2) * dth2).^2;

T = 0.5 * (M * dx.^2 + m1 * v1_sq + m2 * v2_sq);

%---Solve for Potential Energy
V = - (m1 * g * l1 * cos(theta1) ) - (m2 * g * l2 * cos(theta2));

%--Solve for Langrangian of the System
L = T - V;

%%%%%%%%%%%%%%%%%%%%%%%

%---Euler-Lagrange Method for Equations of Motion

dL_x = diff(L, x);
dL_dx = diff(L, dx);
term_x = diff(dL_dx) - dL_x - F;

dL_theta1 = diff(L, theta1);
dL_dtheta1 = diff(L, dth1);
term_theta1 = diff(dL_dtheta1) - dL_theta1;

dL_theta2 = diff(L, theta2);
dL_dtheta2 = diff(L, dth2);
term_theta2 = diff(dL_dtheta2) - dL_theta2;

ddx = diff(dx);
ddth1 = diff(dth1);
ddth2 = diff(dth2);

%-------Solve for DDx, DDtheta1, DDtheta2----------
%---Make Equations Symbolic for solver---
syms Dx Dth1 Dth2 DDx DDtheta1 DDtheta2 x_a th1 th2
term_x = subs (term_x, {x, theta1, theta2, dx, dth1, dth2, ddx, ddth1, ddth2},...
    {x_a, th1, th2, Dx, Dth1, Dth2, DDx, DDtheta1, DDtheta2});

term_theta1 = subs (term_theta1, {x, theta1, theta2, dx, dth1, dth2, ddx, ddth1, ddth2},...
    {x_a, th1, th2, Dx, Dth1, Dth2, DDx, DDtheta1, DDtheta2});

term_theta2 = subs (term_theta2, {x, theta1, theta2, dx, dth1, dth2, ddx, ddth1, ddth2},...
    {x_a, th1, th2, Dx, Dth1, Dth2, DDx, DDtheta1, DDtheta2});

eqns = [term_x; term_theta1; term_theta2];

vars = [DDx; DDtheta1; DDtheta2];


S = solve(eqns, vars);
disp('Equations of Motion')
DDx = S.DDx    %x DD = doubledot
DDtheta1 = S.DDtheta1
DDtheta2 = S.DDtheta2
%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Solve for Matrix A and B using Jacobian---
X = [x_a; th1; th2; Dx; Dth1; Dth2];
Xdot = [Dx; Dth1; Dth2; DDx; DDtheta1; DDtheta2];

U = F;  %Input

J_U = jacobian(Xdot, U);
J_X = jacobian(Xdot, X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------Part B: Plug in Equilibrium Points
% x_a, th1, th2, Dx, Dth1, Dth2
disp('Linearized System about the equilibrium Point')
A = subs (J_X, {sym('x_a'),sym('th1'), sym('th2'), sym('Dx'),sym('Dth1'), sym('Dth2')}, ...
    {0,0,0,0,0,0})

% A =
%  
% [0,                    0,                    0, 1, 0, 0]
% [0,                    0,                    0, 0, 1, 0]
% [0,                    0,                    0, 0, 0, 1]
% [0,             (g*m1)/M,             (g*m2)/M, 0, 0, 0]
% [0, -(M*g + g*m1)/(M*l1),       -(g*m2)/(M*l1), 0, 0, 0]
% [0,       -(g*m1)/(M*l2), -(M*g + g*m2)/(M*l2), 0, 0, 0]

B = subs (J_U,{sym('x_a'),sym('th1'), sym('th2'), sym('Dx'),sym('Dth1'), sym('Dth2')}, ...
    {0,0,0,0,0,0} )
% B =
%  
%         0
%         0
%         0
%       1/M
% -1/(M*l1)
% -1/(M*l2)

%----------------------------------
%----------Part C----------
% Create Controllability Matrix/conditions for
%       M, m1, m2, l1, l2

% C = [C = B AB A2B A3B A4B A5B]

C = B;
temp = B;

n = size(A)-1;

for i = 1:n
%     temp = (A.^i)*temp;
    temp = A*temp;
    C = horzcat(C, temp);
end

disp('Controllability Matrix (C) is ')
C
C_det = simplify (det(C))

% Proves that will be singular when l1=l2
singular_test = solve(C_det==0, [l1]); 

%----------------------------------
%----------Part D--------------------
% Plug in M, m1, m2, l1, l2 values into C matrix
C = subs(C, {M, m1, m2, l1, l2, g},{1000, 100, 100, 20, 10, 10})

%Check for controllability
C_det = sprintf('%10e',det(C)) %'1.562500e-24'

if rank(C) == size(A,1)
    disp('System is controllable!')
else
    disp('System is NOT controllable!!!')
end

%Obtain LQR controller
syms F X_1 X_2 X_3 X_1d X_2d X_3d;
tspan = 0:0.1:200;
% 10 degrees = pi/18
s0 = [0; pi/18; pi/18; 0; 0; 0]; %1 meter push to serve as step input
X = [X_1; X_2; X_3; X_1d; X_2d; X_3d;];

A = subs(A, {M, m1, m2, l1, l2, g},{1000, 100, 100, 20, 10, 10});
B = subs(B, {M, m1, m2, l1, l2, g},{1000, 100, 100, 20, 10, 10});

%Convert sym matrixes to num
A = double (A);
B = double (B);

R = .1; %Bad/initial R
Q = diag([1 1 1 1 1 1]);    %%Bad/initial Q

%Xd = A*X + B*F;
[K,S,e] = lqr(A,B,Q,R);

%-----------------------------------------------------------
%------Simulate resulting responses to I.C.----------

%----No input (without forces acting on it)----
[t,state_history] = ode45(@(t,state)no_input_lin_model(state,t,A),tspan,s0);

figure('Name', 'Linearized Model No Forces');
subplot(2,2,1);
plot(t,state_history(:,1),'k')
grid on;
ylabel('x position of the cart')
xlabel('time in s')

subplot(2,2,2);
plot(t,state_history(:,2),'r')
grid on;

subplot(2,2,3);
plot(t,state_history(:,3),'b')
grid on;

%------LQR Control Signal----
[t,state_history] = ode45(@(t,state)linearized_model(state,t,A,B,K),tspan,s0);

figure('Name', 'Linearized Model with Force (Bad Q&R Values)');
subplot(2,2,1);
plot(t,state_history(:,1),'k')
grid on;
ylabel('x position of the cart')
xlabel('time in s')

subplot(2,2,2);
plot(t,state_history(:,2),'r')
grid on;

subplot(2,2,3);
plot(t,state_history(:,3),'b')
grid on;

%----------------------------
% Original Non-linear system
[t,state_history] = ode45(@(t,state)og_model(state,t,K),tspan,s0);

figure('Name', 'Non-Linearized Model with Force (Bad Q&R Values)');
subplot(2,2,1);
plot(t,state_history(:,1),'k')
grid on;
ylabel('x position of the cart')
xlabel('time in s')

subplot(2,2,2);
plot(t,state_history(:,2),'r')
grid on;

subplot(2,2,3);
plot(t,state_history(:,3),'b')
grid on;

%---------------------------------------
%------Adjust parameters of LQR cost
%---------------------------------
R = .00001;
Q = diag([10 1000 1000 1000 1000 10000]);

[K,S,e] = lqr(A,B,Q,R);
disp('LQR Controller control Matrix (K) is ')
K
%K = [0.1000   -1.0473   -3.2167    1.0532   -0.6226   -0.0986]

[t,state_history] = ode45(@(t,state)linearized_model(state,t,A,B,K),tspan,s0);

figure('Name', 'Linearized Model with Force (Good Q&R Values)');
subplot(2,2,1);
plot(t,state_history(:,1),'k')
grid on;
ylabel('x position of the cart')
xlabel('time in s')

subplot(2,2,2);
plot(t,state_history(:,2),'r')
grid on;

subplot(2,2,3);
plot(t,state_history(:,3),'b')
grid on;

%----------------------------
% Original Non-linear system
[t,state_history] = ode45(@(t,state)og_model(state,t,K),tspan,s0);

figure('Name', 'Non-Linearized Model with Force (Good Q&R Values)');
subplot(2,2,1);
plot(t,state_history(:,1),'k')
grid on;
ylabel('x position of the cart')
xlabel('time in s')

subplot(2,2,2);
plot(t,state_history(:,2),'r')
grid on;

subplot(2,2,3);
plot(t,state_history(:,3),'b')
grid on;

%Lyapunov's indirect method to certify Stability of closed loop system
disp('Eigenvalues using Lyapunov indirect method')
Lyp = eig(A-B*K)
stability_check = 1;
% for i = 1:6
for i = 1:size(A,1)
    if real(Lyp(i)) >= 0
            stability_check = 0;
        break;
    end
end
if stability_check == 1
    disp("Eigen values all have negative, real parts, so it is stable locally");
end

%----------------------------------
%----------Second Component--------------------
%---Part E----------------------
%Output Vectors
disp('C matrices')
V1 = [1 0 0 0 0 0]             %For output (x(t))
V2 = [0 1 0 0 0 0;           %For output (t1(t),th2(t)) 
      0 0 1 0 0 0]
V3 = [1 0 0 0 0 0;           %For output (x(t),th2(t))    
      0 0 1 0 0 0]
V4 = [1 0 0 0 0 0;           %For output (x(t),th1(t),th2(t))
      0 1 0 0 0 0;
      0 0 1 0 0 0]

%Check Observability (rank = n)

obs_test1 = [V1;V1*A;V1*A^2;V1*A^3;V1*A^4;V1*A^5;];
obs_test2 = [V2;V2*A;V2*A^2;V2*A^3;V2*A^4;V2*A^5;];
obs_test3 = [V3;V3*A;V3*A^2;V3*A^3;V3*A^4;V3*A^5;];
obs_test4 = [V4;V4*A;V4*A^2;V4*A^3;V4*A^4;V4*A^5;];

% Calculate ranks for Observability
% disp('Ranks')
obs_res1 = rank(obs_test1);
obs_res2 = rank(obs_test2);
obs_res3 = rank(obs_test3);
obs_res4 = rank(obs_test4);


if obs_res1 == size(A,1)
    disp("System is observable for x(t)");  %True
else
    disp("System is NOT observable for x(t)");
end
if obs_res2 == size(A,1)
    disp("System is observable for th1(t),th2(t))");
else
    disp("System is NOT observable for th1(t),th2(t))"); %True
end
if obs_res3 == size(A,1)
    disp("System is observable for x(t),th2(t)");   %True
else
    disp("System is NOT observable for x(t),th2(t)");
end
if obs_res4 == size(A,1)
    disp("System is observable for x(t,th1(t),th2(t))");    %True
else
    disp("System is NOT observable for x(t,th1(t),th2(t))");
end

%----------------------------------
%--------Part F-------
%Luenberger Observer for x(t), (x(t),th2(t)), (x(t), th1(t),th2(t))

D=[0]; %Set to 0
p=[-2 -3 -4 -5 -6 -7 ]; % Needed values for p to meet place rank requirement
x_0= [ 0 ; 0 ; 0.1 ; 0.1 ;0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0];

%---For output (x(t))---
figure('Name', 'Luenberger Observer for x(t)')
C = V1;

% Pole Placement method
L=place(A',C',p);
L=L';

%New State Matrices with observer
Ac=[(A-B*K) (B*K); zeros(size(A)) (A-L*C)];
Bc=[ B ;zeros(size(B))];
Cc= [C zeros(size(C))];
Dc=D;

sys=ss(Ac,Bc,Cc,Dc);
step(sys);
initial(sys,x_0);
grid on

%---For output (x(t),th2(t)) ---
figure('Name', 'Luenberger Observer for x(t),th2(t)')
C = V3;
L=place(A',C',p);
L=L';

Ac=[(A-B*K) (B*K); zeros(size(A)) (A-L*C)];
Cc= [C zeros(size(C))];

sys=ss(Ac,Bc,Cc,Dc);
step(sys);
initial(sys,x_0)
grid on

%---For output (x(t),t1(t),t2(t))---
figure('Name', 'Luenberger Observer for x(t),th1(t),th2(t)')
C = V4;

L=place(A',C',p);
L=L';

Ac=[(A-B*K) (B*K); zeros(size(A)) (A-L*C)];
Cc= [C zeros(size(C))];

sys=ss(Ac,Bc,Cc,Dc);
step(sys);
initial(sys,x_0)
grid on

%----------------------------------
%----------Part G--------
%--Design LQG controller to original nonlinear system
% X_dot = AX + BU + BW = CX + V

% Noise Variables
Bd = 0.001.*eye(6);             %input disturbance covariance
Bn1 = 0.01;                      %output measurement noise
Bn3 = 0.01*eye(2);
Bn4 = 0.01*eye(3);

% Obtaining Luenberger observers for use in LQG
[L1,P,E] = lqe(A,Bd,V1,Bd,Bn1);
[L3,P,E] = lqe(A,Bd,V3,Bd,Bn3);
[L4,P,E] = lqe(A,Bd,V4,Bd,Bn4);

%x(t) is best "smallest" output vector

%LQG Section for both linear and non-linear using x(t) output vector

%Linearized system to make LQG using 
C = V1; %x(t) was best "smallest" output vector
sys_lqr = ss(A,[B B],C,[zeros(1,1) zeros(1,1)]);
vd = 0.1 * eye(1);
vn = 1;
sen = [1];
known = [1];

%Obtain l Matrix using Kalman filter
[~,L,~] = kalman(sys_lqr,vd,vn,[],sen,known)
states = {'x','theta1','theta2','x_dot','theta1_dot','theta2_dot','e_1','e_2','e_3','e_4','e_5','e_6'};
inputs = {'F'};
outputs = {'x'};

% Closed Loop Controller on original nonlinear system
Ac = [A-B*K B*K;zeros(size(A)) A-L*C];
Bc = zeros(12,1);
Cc = [C zeros(size(C))];
sys_cl_lqg = ss(Ac,Bc,Cc,D, 'statename',states,'inputname',inputs,'outputname',outputs);    

% Set initial conditions to x=5; th1=5 degrees, th2=10 degrees
x0 = [ 5 ; 0.1 ; 0.2 ; 0 ;0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0];
t = 0:0.01:100;
F = zeros(size(t));
[Y,~,X] = lsim(sys_cl_lqg,F,t,x0);

% Input vector: u=K * x(t)
u = zeros(size(t));
for i = 1:size(X,1)
   u(i) = K * (X(i,1:6))';
end
Xhat = X(:,1) - X(:,6);

figure ('Name', 'LQG Controller applied to original system')
subplot(2,2,1);
plot(t,Y(:,1),'k')
ylabel('system response of the cart')
xlabel('time in s')
grid on;

subplot(2,2,2);
plot(t,Xhat, 'b')
ylabel('xHat of the cart')
xlabel('time in s')
grid on;

subplot(2,2,3);
plot(t,X(:,1),'r')
ylabel('x position of the cart')
xlabel('time in s')
grid on;


function sdot = linearized_model(s,t,A,B,K)
%% Linearized Model
sdot = (A-B*K)*s;
end


function sdot = og_model(s,t,K)
%% Non-Linearized Model
F = -K*s;
sdot = [s(4); s(5); s(6);
    (F + 2000*sin(s(2))*s(5)^2 + 500.0*sin(2.0*s(2)) + 1000*sin(s(3))*s(6)^2 + 500.0*sin(2.0*s(3)))/(100*sin(s(2))^2 + 100*sin(s(3))^2 + 1000);
    -(10*(100*sin(s(3)))^2 + 1100)*sin(s(2)) + (F + 2000*sin(s(2))*s(5)^2 + 1000*sin(s(3))*s(6)^2)*cos(s(2)) - 250.0*sin(s(2) - 2.0*s(3)) + 250.0*sin(s(2) + 2.0*s(3))/(20*(100*sin(s(2))^2 + 100*sin(s(3))^2 + 1000));
    -(10*(100*sin(s(2))^2 + 1100)*sin(s(3)) + (F + 2000*sin(s(2))*s(5)^2 + 1000*sin(s(3))*s(6)^2)*cos(s(3)) + 250.0*sin(2.0*s(2) - s(3)) + 250.0*sin(2.0*s(2) + s(3)))/(10*(100*sin(s(2))^2 + 100*sin(s(3))^2 + 1000));     
];
end 


function sdot = no_input_lin_model(s,t,A)
%% Model without forces acting on it
sdot = (A)*s;
end

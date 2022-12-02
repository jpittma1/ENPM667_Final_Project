%% ENPM 667 Fall 2022
%% Final Project
%% Jerry Pittman (117077120
%% Bob Reiter (UID)

clear; 
clc;

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


S = solve(eqns, vars)

DDx = S.DDx;    %x doubledot
DDtheta1 = S.DDtheta1;
DDtheta2 = S.DDtheta2;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Solve for Matrix A and B using Jacobian---
X = [x_a; th1; th2; Dx; Dth1; Dth2];
Xdot = [Dx; Dth1; Dth2; DDx; DDtheta1; DDtheta2];

U = [F];  %Input

J_U = jacobian(Xdot, U)
J_X = jacobian(Xdot, X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------Part B: Plug in Equilibrium Points
% x_a, th1, th2, Dx, Dth1, Dth2
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

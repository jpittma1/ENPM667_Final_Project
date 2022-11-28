clear;
clc;



syms F X_1 X_2 X_3 X_1d X_2d X_3d;
tspan = 0:0.1:200;
s0 = [0; pi/18; pi/18; 0; 0; 0]; %1 meter push to serve as step input
X = [X_1; X_2; X_3; X_1d; X_2d; X_3d;];

A = [0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    0 1.00000000000000 1.00000000000000 0 0 0;
    0 -11/20 -0.0500000000000000 0 0 0;
    0 -0.100000000000000 -11/10 0 0 0;
    ];

B = [0; 0; 0; 1/1000; -1/20000; -1/10000;];
R = .00001;
Q = diag([10 1000 1000 1000 1000 10000]);

[K,S,e] = lqr(A,B,Q,R);


%Xd = A*X + B*F;


[t,state_history] = ode45(@(t,state)no_input_lin_model(state,t,A),tspan,s0);

figure;
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

[t,state_history] = ode45(@(t,state)linearized_model(state,t,A,B,K),tspan,s0);

figure;
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

A = [0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    0 10.000000000 1.00000000000000 0 0 0;
    0 -101/20 -0.0500000000000000 0 0 0;
    0 -.1000000000000 -101/10 0 0 0;
    ];

B = [0; 0; 0; 1/1000; -1/20000; -1/10000;];

tspan = 0:0.1:200;
[t,state_history] = ode45(@(t,state)linearized_model(state,t,A,B,K),tspan,s0);

figure;
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

tspan = 0:0.1:200;


A = [0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    0 1.00000000000000 1.00000000000000 0 0 0;
    0 -11/20 -0.0500000000000000 0 0 0;
    0 -0.100000000000000 -11/10 0 0 0;
    ];

B = [0; 0; 0; 1/1000; -1/20000; -1/10000;];




% Lyupanov's Indirect Method

Lyp = eig(A-B*K);
stability_check = 1;
for i = 1:6
    if real(Lyp(i)) >= 0
            stability_check = 0;
        break;
    end
end
if stability_check == 1
    disp("Eigen values all have negative, real parts, so it is stable locally");
end



% Part 2



V1 = [1 0 0 0 0 0];             %For output (x(t))
V2 = [0 1 0 0 0 0;           %For output (t1(t),t2(t)) 
      0 0 1 0 0 0];
V3 = [1 0 0 0 0 0;           %For output (x(t),t2(t))    
      0 0 1 0 0 0];
V4 = [1 0 0 0 0 0;           %For output (x(t),t1(t),t2(t))
      0 1 0 0 0 0;
      0 0 1 0 0 0];
  
% Observability Confirmation:
% Is rank = n?

ob_m1 = [V1;V1*A;V1*A^2;V1*A^3;V1*A^4;V1*A^5;];
ob_m2 = [V2;V2*A;V2*A^2;V2*A^3;V2*A^4;V2*A^5;];
ob_m3 = [V3;V3*A;V3*A^2;V3*A^3;V3*A^4;V3*A^5;];
ob_m4 = [V4;V4*A;V4*A^2;V4*A^3;V4*A^4;V4*A^5;];

obs_res1 = rank(ob_m1);
obs_res2 = rank(ob_m1);
obs_res3 = rank(ob_m1);
obs_res4 = rank(ob_m1);


if obs_res1 == 6
    disp("System is observable for x(t)");
end
if obs_res2 == 6
else
    disp("System is not observable for th1(t),th2(t))");
end
if obs_res3 == 6
    disp("System is observable for x(t),th2(t)");
end
if obs_res4 == 6
    disp("System is observable for x(t,th1(t),th2(t))");
end

% First vector:
%Setup new figure
figure;

% Step response:
C=[1 0 0 0 0 0 ]; %x(t)
D=[0]; %Set to 0
p=[-2 -3 -4 -5 -6 -7 ]; % Needed values for p to meet place rank requirement
x_0= [ 0 ; 0 ; 0.1 ; 0.1 ;0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0];

% Placement method
L=place(A',C',p);
L=L';

% Define new state matrices
Ac=[(A-B*K) (B*K); zeros(size(A)) (A-L*C)];
Bc=[ B ;zeros(size(B))];
Cc= [C zeros(size(C))];
Dc=[0];

sys_amp_1=ss(Ac,Bc,Cc,Dc);
step(sys_amp_1);

% output section:

initial(sys_amp_1,x_0);

figure;

%Second vector:
C=[1 0 0 0 0 0 ; 0 0 0 0 1 0];
D=[0];


p=[-2 -3 -4 -5 -6 -7 ];
L=place(A',C',p);
L=L';
Ac=[(A-B*K) (B*K); zeros(size(A)) (A-L*C)];
Bc=[ B ;zeros(size(B))];
Cc= [C zeros(size(C))];
Dc=[0];
sys=ss(Ac,Bc,Cc,Dc);
step(sys);

%output section:

initial(sys,x_0)
grid on



figure;
%Last Vector:

C=[1 0 0 0 0 0 ; 0 0 1 0 0 0; 0 0 0 0 1 0]; % Matrix representation of vector
D=[0];

p=[-2 -3 -4 -5 -6 -7 ]
L=place(A',C',p)
L=L'
Ac=[(A-B*K) (B*K); zeros(size(A)) (A-L*C)];
Bc=[ B ;zeros(size(B))];
Cc= [C zeros(size(C))];
Dc=[0];
sys=ss(Ac,Bc,Cc,Dc);
step(sys);


initial(sys,x_0);


% Noise Variables
Bd = 0.001.*eye(6);             %input disturbance covarianve
Bn1 = 0.01;                      %output measurement noise
Bn3 = 0.01*eye(2);
Bn4 = 0.01*eye(3);

% Obtaining Luenberger observers for use in LQG
[L1,P,E] = lqe(A,Bd,V1,Bd,Bn1);
[L3,P,E] = lqe(A,Bd,V3,Bd,Bn3);
[L4,P,E] = lqe(A,Bd,V4,Bd,Bn4);
%LQG Section for both linear and non-linear
A_c = A-L1*V1;
B_c = [B L1];
C_c = eye(6);
D_c = 0*[B L1];

nx = 3;
ny = 1;
R_n = [1 0;0 2];
QXU = eye(14);
QWV = eye(8);
QI = eye(ny);


D = 0;
C = [1 0 0 0 0 0];
sys_lqr = ss(A,[B B],C,[zeros(1,1) zeros(1,1)]);
vd = 0.1 * eye(1);
vn = 1;
sen = [1];
known = [1];
[~,L,~] = kalman(sys_lqr,vd,vn,[],sen,known)
states = {'x','theta1','theta2','x_dot','theta1_dot','theta2_dot','e_1','e_2','e_3','e_4','e_5','e_6'};
inputs = {'F'};
outputs = {'x'};

Ac = [A-B*K B*K;zeros(size(A)) A-L*C];
Bc = zeros(12,1);
Cc = [C zeros(size(C))];
sys_cl_lqg = ss(Ac,Bc,Cc,D, 'statename',states,'inputname',inputs,'outputname',outputs);    
x0 = [ 5 ; 0.1 ; 0.2 ; 0 ;0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0];
t = 0:0.01:100;
F = zeros(size(t));
[Y,~,X] = lsim(sys_cl_lqg,F,t,x0);
figure
plot(t,Y(:,1),'b');

u = zeros(size(t));
for i = 1:size(X,1)
   u(i) = K * (X(i,1:6))';
end
Xhat = X(:,1) - X(:,6);
figure
subplot(3,1,1), plot(t,Xhat), hold on, plot(t,X(:,1),'r')




function sdot = linearized_model(s,t,A,B,K)
%% Linearized Model
sdot = (A-B*K)*s;
end


function sdot = og_model(s,t,K)
%% Linearized Model
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


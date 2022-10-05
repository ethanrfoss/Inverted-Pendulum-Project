%% Inverted Pendulum Controls

%% System:
% Environment:
S.g = 9.81;

% System Characteristics:
encoderMass = 142; %g
sliderMass = 200; %g
beltMass = 20; %g
encoderPlateMass = 15; %g
gantryPlateMass = 30; %g
pendDensity = 2700;
pendDiam = .009525; %m
pendLength = .6096;

S.mc = (encoderMass+sliderMass+beltMass+encoderPlateMass+gantryPlateMass)/1000; % Cart mass[kg]
S.m = pi*(pendDiam/2)^2*pendLength*pendDensity; % Rod Mass[kg]
S.R = pendLength; % Rod Length[m]
S.l = S.R/2;
S.st = .3358;
S.k = .0046;
S.I = pi/2*2700*(7/1000*(25/1000)^4-(6/1000)^4+10/1000*((42/1000)^4-(6/1000)^4));
S.r = 21/1000;

S.B = [0;0;2*S.st/S.r;0];
%% Observability, Controllability
syms a32 a34 a42;
A = [0 0 1 0;0 0 0 0;0 -a32 0 -a34;0 a42 0 0];
B = [0;0;1;0];
C = [1 0 0 0;0 1 0 0];

Con = [B A*B A^2*B A^3*B];
Obs = [C; C*A; C*A^2; C*A^3];
%% MPC:
Q=diag([0 0 0 0]); R=10000; QT=diag([1 5000 2 800000]); T = 3;
[u,x,t] = IPMPC(T,Q,R,QT);
x(:,end)
%invertedPendulumMPCTester(t,u,x(:,1))
% S.Q=diag([0 10 0 0]);S.R=0;S.QT=diag([2000 8000 0 100]);
% S.Q=diag([0 0 0 0]);S.R=0;S.QT=diag([20 80 0 10]);
% S.Q=diag([0 0 0 0]); S.R=0; S.QT=diag([20 100 0 5]);
% S.Q=diag([0 0 0 0]); S.R=1;S.QT=diag([20 100 0 20]);
% Q=diag([0 10 0 0]);R=50;QT=diag([20 100 0 30]);
% Q=diag([0 10 0 0]); R=50; QT=diag([0 100 0 7000]);
% Q=diag([0 10 0 0]);R=50;QT=diag([1 200 4 6000]);
% Q=diag([0 10 0 0]); R=50; QT=diag([1 500 3 7000]);
% Q=diag([0 10 0 0]); R=50; QT=diag([1 550 2 7300]);
% Q=diag([0 10 0 0]); R=50; QT=diag([1 520 2 7300]);
% Q=diag([10 0 0 0]); R=1000; QT=diag([1 520 2 5000]); T = 3;
% Q=diag([0 0 0 0]); R=10000; QT=diag([1 5000 2 800000]); T = 3; k = 300

%% LQR infinite horizon:
x0 = [0;-pi;0;0];
Qk = diag([1000;1000;0;100]); Rk = 1000;
Ql = diag([0.000001,0.00001,0,0.0000000001]); Rl = diag([.0052,1.0996*10^-4]);
%invertedPendulumLQRTester(x0,Q,R);
%invertedPendulumLQGTester(x0,Qk,Rk,Ql,Rl);
[Kk,K,xk] = invertedPendulumMPCLQRTester(x0,Qk,Rk,u,t,x);


function invertedPendulumLQRTester(x0,Q,R)
%% Double Pendulum Simulation:
%% Environment:
S.g = 9.81;

%% System Characteristics:
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

%% Gain Matrix:
X = icare(computeA(S,[0;0;0;0;0;0]),S.B,Q,R,0,computeE(S,[0;0;0;0]));
K = -R^-1*S.B'*X*computeE(S,[0;0;0;0])

%% Initial State:
x(:,1) = x0(1:4);

%% Shapes for plotting:
warning off;
bodyI = polyshape([.5*S.R .5*S.R -.5*S.R -.5*S.R],[.25*S.R -.25*S.R -.25*S.R .25*S.R]);
thickness = S.R/32;
pendI = polyshape([thickness thickness -thickness -thickness],[0 S.R S.R 0]);
body = bodyI;
pend = pendI;

%% Loop Signals:
uCom = 0;
tCom = 0;
et = 0;
ex = 0;

%% Loop:

h = figure; hold on; axis equal; grid on; tit = title({'Double Pendulum System',sprintf('t = 0')});
xlim([min(x(1,1))-2*S.R-.2 max(x(1,1))+2*S.R+.2]); ylim([-2*S.R-.2, 2*S.R+.2]);
plot([-100-2*S.R-.2 100+2*S.R+.2],[0 0]);
pb = plot(bodyI);
pp = plot(pendI);

dt = .01;
t = [0:dt:40];

for i = 1:length(t)
    
    uCom = K*x(:,i);
    
    fprintf("| Motor Com: %.2f | Angle: %.2f  | Position: %.2f |\n",uCom,x(2,i),x(1,i));
    
    x(:,i+1) = RK4(@computeE,@computeN,S,x(:,i),uCom,dt);
    
    body.Vertices = bodyI.Vertices +[x(1,i+1) 0];
    pend.Vertices = rotate(pendI,x(2,i+1)*180/pi).Vertices + [x(1,i+1) 0];
    
    xlim([min(x(1,i+1))-2*S.R-.2, max(x(1,i+1))+2*S.R+.2]);
    set(pb,'Shape',body);
    set(pp,'Shape',pend);
    drawnow;
    
    set(tit,'String',{'Double Pendulum System',sprintf('t = %.2f',t(i))})

end

end

function E = computeE(S,x)

E = [eye(2,2), zeros(2,2); zeros(2,2), [(S.mc+S.m+2*S.I/S.r^2), -S.m*S.l*cos(x(2));-S.m*S.l*cos(x(2)),4/3*S.m*S.l^2]];

end

function A = computeA(S,x)

A = [zeros(2,2), eye(2,2);[0,-S.m*S.l*(x(6)*sin(x(2))+x(4)^2*cos(x(2)));0,S.m*S.l*(S.g*cos(x(2))-x(5)*sin(x(2)))],[-2*S.k/S.r^2,-2*S.m*S.l*x(4)*sin(x(2));0,0]];

end

function N = computeN(S,x,u)

if u < .07 && u > -.07
    u = 0;
end
if u < -1
    u = -1;
end
if u > 1
    u = 1;
end

N = [x(3);x(4);2*(S.st*S.r*u-S.k*x(3))/S.r^2 - S.m*S.l*sin(x(2))*x(4)^2;S.m*S.g*S.l*sin(x(2))];

end

function xNext = RK4(computeE,computeN,S,x,u,dt)

f1 = computeE(S,x)\computeN(S,x,u);
f2 = computeE(S,x+dt*f1/2)\computeN(S,x+dt*f1/2,u);
f3 = computeE(S,x+dt*f2/2)\computeN(S,x+dt*f2/2,u);
f4 = computeE(S,x+dt*f3)\computeN(S,x+dt*f3,u);

xNext = x + dt*(f1/6+(f2+f3)/3+f4/6);

end

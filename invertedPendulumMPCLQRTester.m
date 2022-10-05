function [Kk,K,x] = invertedPendulumMPCLQRTester(x0,Q,R,uk,tk,xk)
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

%% Finite Horizon LQR:
h = tk(2)-tk(1);
for i = length(tk):-1:1
    Kk{i} = -R^-1*S.B'*X*computeE(S,xk(:,i));
    
    f1=RHSX(xk(:,i),S,X,Q,R);
    f2=RHSX(xk(:,i),S,X-h*f1/2,Q,R);
    f3=RHSX(xk(:,i),S,X-h*f2/2,Q,R);
    f4=RHSX(xk(:,i),S,X-h*f3,Q,R);
    
    X=X-h*(f1/6+(f2+f3)/3+f4/6);
end

%% Initial State:
x(:,1) = x0(1:4);
xm = [x0(1:2); zeros(2,1); zeros(2,1)];

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
t = [0:dt:10];

for i = 1:length(t)
    
    if mod(t(i),.02) == 0
        if t(i) > tk(end)
            uCom = K*xm(1:4);
        else
            n = find(min(abs(t(i)-tk)) == t(i)-tk);
            uCom = uk(n) + Kk{n}*(xm(1:4)-xk(1:4,n));
        end
    end
    
    if mod(t(i),.02) == 0
        xm = updateState2(xm,x(:,i),.02);
    end
    
    
    %fprintf("| Motor Com: %.2f | Angle: %.2f  | Position: %.2f |\n",uCom,x(2,i),x(1,i));
    fprintf("| Motor Com: %.2f | t: %.2f | Pos: %.2f %.2f | Ang: %.2f %.2f | Pos Rate: %.2f %.2f | Ang Rate: %.2f %.2f |\n",uCom,t(i),x(1,i),xm(1),x(2,i),xm(2),x(3,i),xm(3),x(4,i),xm(4));
    x(:,i+1) = RK4(@computeE,@computeN,S,x(:,i),uCom,dt);
    
    body.Vertices = bodyI.Vertices +[x(1,i+1) 0];
    pend.Vertices = rotate(pendI,x(2,i+1)*180/pi).Vertices + [x(1,i+1) 0];
    
    xlim([min(x(1,i+1))-2*S.R-.2, max(x(1,i+1))+2*S.R+.2]);
    set(pb,'Shape',body);
    set(pp,'Shape',pend);
    drawnow;
    
    set(tit,'String',{'Double Pendulum System',sprintf('t = %.2f',t(i))})
    
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if i == 1
%         imwrite(imind,cm,[cd '\SavedGifs\' 'MPCSim' '.gif'],'gif','DelayTime',0, 'Loopcount',inf); 
%     else
%         imwrite(imind,cm,[cd '\SavedGifs\' 'MPCSim' '.gif'],'gif','DelayTime',0,'WriteMode','append'); 
%     end

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

function xm = updateState(xm,x,dt)

encoder = [-1200:1200]*2*pi/1200;
pos = [-1:2*pi/1200*21/1000:1];

newAng = encoder(find(abs(x(2)-encoder) == min(abs(x(2)-encoder))));
newPos = pos(find(abs(x(1)-pos) == min(abs(x(1)-pos))));

xm(4) = (newAng-xm(2))/dt;
xm(3) = (newPos-xm(1))/dt;
xm(2) = newAng;
xm(1) = newPos;

end


function xm = updateState2(xm,x,dt)

encoder = [-1200:1200]*2*pi/1200;
pos = [-1:2*pi/1200*21/1000:1];

newAng = encoder(find(abs(x(2)-encoder) == min(abs(x(2)-encoder))));
newPos = pos(find(abs(x(1)-pos) == min(abs(x(1)-pos))));

if newAng == xm(2)
    xm(6) = xm(6) + dt;
else
    xm(4) = (newAng -xm(2))/(xm(6)+dt);
    xm(6) = 0;
end

if newPos == xm(1)
    xm(5) = xm(5) + dt;
else
    xm(3) = (newPos -xm(1))/(xm(5)+dt);
    xm(5) = 0;
end

xm(2) = newAng;
xm(1) = newPos;

end

function R=RHSX(x,S,X,Q,R)
E = computeE(S,x);
A = computeA(S,x);
R = -(E')^-1*A'*X-X*A*E^-1+X*S.B*R^-1*S.B'*X-(E')^-1*Q*E^-1;
end
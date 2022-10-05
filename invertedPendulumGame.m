function invertedPendulumGame
%% Double Pendulum Simulation:
%% Environment:
S.g = 9.81;

%% System Characteristics:
S.mc = .189*4+0.141748; % Cart mass[kg]
S.m = .4692; % Rod Mass[kg]
S.R = .6096; % Rod Length[m]
S.l = S.R/2;

%% Initial State:
x(:,1) = [0;... %x
          pi;... %theta
          0;... %xdot
          0]; %thetadot

%% Shapes for plotting:
warning off;
bodyI = polyshape([.5*S.R .5*S.R -.5*S.R -.5*S.R],[.25*S.R -.25*S.R -.25*S.R .25*S.R]);
thickness = S.R/32;
pendI = polyshape([thickness thickness -thickness -thickness],[0 S.R S.R 0]);
body = bodyI;
pend = pendI;

%% Loop:

h = figure; hold on; axis equal; grid on; tit = title({'Double Pendulum System',sprintf('t = 0')});
xlim([min(x(1,1))-2*S.R-.2 max(x(1,1))+2*S.R+.2]); ylim([-2*S.R-.2, 2*S.R+.2]);
plot([-100-2*S.R-.2 100+2*S.R+.2],[0 0]);
pb = plot(bodyI);
pp = plot(pendI);

dt = .025;
run = true;
t = 0;
i = 1;

while run
    
    keysdown = keys();
    F = 0;
    [n,~] = size(keysdown);
    for j = 1:n
        if isequal('Left',keysdown{j})
            F = -100;
        elseif isequal('Right',keysdown{j})
            F = 100;
        elseif isequal('Enter',keysdown{j})
            run = false;
        end
    end
    
    x(:,i+1) = RK4(@computeE,@computeN,S,x(:,i),F,dt);
    t = t + dt;
    
    body.Vertices = bodyI.Vertices +[x(1,i+1) 0];
    pend.Vertices = rotate(pendI,x(2,i+1)*180/pi).Vertices + [x(1,i+1) 0];
    
    xlim([min(x(1,i+1))-2*S.R-.2, max(x(1,i+1))+2*S.R+.2]);
    set(pb,'Shape',body);
    set(pp,'Shape',pend);
    drawnow;
    
    set(tit,'String',{'Double Pendulum System',sprintf('t = %.2f',t)});
    
    i = i + 1;

end

end
function E = computeE(S,x)

E = [eye(2,2), zeros(2,2); zeros(2,2), [(S.mc+S.m), -S.m*S.l*cos(x(2));-S.m*S.l*cos(x(2)),4/3*S.m*S.l^2]];

end

function N = computeN(S,x,F)

N = [x(3);x(4);F-S.m*S.l*sin(x(2))*x(4)^2;S.m*S.g*S.l*sin(x(2))];

end

function xNext = RK4(computeE,computeN,S,x,F,dt)

f1 = computeE(S,x)\computeN(S,x,F);
f2 = computeE(S,x+dt*f1/2)\computeN(S,x+dt*f1/2,F);
f3 = computeE(S,x+dt*f2/2)\computeN(S,x+dt*f2/2,F);
f4 = computeE(S,x+dt*f3)\computeN(S,x+dt*f3,F);

xNext = x + dt*(f1/6+(f2+f3)/3+f4/6);

end

function keysdown = keys()
NET.addAssembly('PresentationCore');
akey = System.Windows.Input.Key.A;  %use any key to get the enum type
keys = System.Enum.GetValues(akey.GetType);  %get all members of enumeration
keynames = cell(System.Enum.GetNames(akey.GetType))';
iskeyvalid = true(keys.Length, 1);
iskeydown = false(keys.Length, 1);
keysdown = [];
for keyidx = 1:keys.Length
   try
       iskeydown(keyidx) = System.Windows.Input.Keyboard.IsKeyDown(keys(keyidx));
   catch
       iskeyvalid(keyidx) = false;
   end
end
keysdown = keynames(iskeydown);
end


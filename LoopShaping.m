%% Loop Shaping Control Design With Simple Model

%% System Paramters

encoderMass = 142; %g
sliderMass = 200; %g
beltMass = 20; %g
encoderPlateMass = 15; %g
gantryPlateMass = 30; %g

pendLength = .6096; %m
pendDiam = .009525; %m
pendDensity = 2700; %kg/m^3

st = .3358;
k = .0046;

I = pi/2*2700*(7/1000*(25/1000)^4-(6/1000)^4+10/1000*((42/1000)^4-(6/1000)^4));
M = (encoderMass+sliderMass+beltMass+encoderPlateMass+gantryPlateMass)/1000;
m = pi*(pendDiam/2)^2*pendLength*pendDensity;
r = 21/1000;
l = pendLength/2;

g = 9.81;

%% System Models

G1 = tf([2*st*r 0],[(((M+m)*r^2+2*I)*4*l/3-m*l*r^2) 8*k*l/3 -g*((M+m)*r^2+2*I) -2*k*g]);% U to Theta

G2 = tf([4*l/3 0 -g],[1 0 0]); % Theta to X

%% Controllers

D1 = 284.29*tf([.18 1],[6.7 1]);
P1 = 1/1.2;

D2= -.017858*tf([1.7 1],[.1 1]);

%% Full System

T1 = minreal(D1*G1/(1+D1*G1));
GD2 = T1*P1*G2*D2;

[num dem] = tfdata(GD2, 'v');
num = tf(num,1);
dem = tf(dem,1);

GT = minreal(num/(dem+num));

%% Discrete Controllers

w1= 11.3;
w2= .472;
D1z=c2d(D1,1/50,['Method','tustin', 'PrewarpFrequency', w1]);
D2z=c2d(D2,1/10,['Method','tustin', 'PrewarpFrequency', w2]);


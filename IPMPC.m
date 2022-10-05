function[uk,xk,t]=IPMPC(T,Q,R,QT,uk)
%STEP0: initialize simulation system & derived parameters
S.h=0.01;
S.N=T/S.h;
S.g = 9.81;
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
t=[0:S.N]*S.h;
alpha=0.1;
S.B=[0;0;2*S.st/S.r;0];
S.Q=Q;
S.R=R;
S.QT=QT;
if nargin<5
    uk=zeros(S.N+1,1);
end
S.x0=[0;-pi;0;0];
xk(1:4,1)=S.x0;
res=0;
for k=0:300
    k
    u=uk(1);
    x=S.x0;
    J=0.25*S.h*(x'*S.Q*x+u'*S.R*u);
    c=.5;
    %STEP 1: march/save state(from t=0−>T), compute cost
    for n=1:S.N
        u=uk(n);
        f1=RHS(x,u,S);
        f2=RHS(x+S.h*f1/2,u,S);
        f3=RHS(x+S.h*f2/2,u,S);
        f4=RHS(x+S.h*f3,u,S);
        x=x+S.h*(f1/6+(f2+f3)/3+f4/6);
        xk(1:4,n+1)=x;
        u=uk(n+1);
        xk(5:6,n)=f1(3:4);
        if n==S.N
            c=.25;
        end
        J=J+c*S.h*(x'*S.Q*x+u'*S.R*u);
    end
    f1=RHS(x,u,S);
    xk(5:6,S.N+1)=f1(3:4);
    E=ComputeE(x,S);
    J=J+0.5*(x'*E'*S.QT*E*x);
    r=S.QT*E*x;
    g(S.N+1,1)=S.B'*r+S.R*uk(S.N+1);
    %STEPS 2&3: march adjoint(from t=T−>0), compute gradient
    for n=S.N:-1:1
        xh=(xk(:,n+1)+xk(:,n))/2;
        f1=RHSa(r,xk(:,n+1),S);
        f2=RHSa(r-S.h*f1/2,xh,S);
        f3=RHSa(r-S.h*f2/2,xh,S);
        f4=RHSa(r-S.h*f3,xk(:,n),S);
        r=r-S.h*(f1/6+(f2+f3)/3+f4/6);
        g(n,1)=S.B'*r+S.R*uk(n);
    end
    res1=res;res=g'*g;
    %STEPS 4&5: update u and repeat
    if(mod(k,4)==0||alpha<1e-4)
        pk=-g;
    else
        pk=-g+pk*res/res1;
    end
    %conjugate gradient
    h = figure(1);clf;    
    subplot(2,1,1);
    plot(t,xk(1,:),'r-',t,xk(2,:),'b-');
    legend('x','\theta1')
    
    subplot(2,1,2);
    plot(t,uk,'r--');
    % Optional Gif
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if k == 0
%         imwrite(imind,cm,[cd '\SavedGifs\' 'MPC' '.gif'],'gif','DelayTime',0, 'Loopcount',inf); 
%     else
%         imwrite(imind,cm,[cd '\SavedGifs\' 'MPC' '.gif'],'gif','DelayTime',0,'WriteMode','append'); 
%     end
    
    [AA,AB,AC,JA,JB,JC]=Bracket(@ComputeJEx211,0,alpha,J,uk,pk,S);%find triplet
    [alpha,J]=Brent(@ComputeJEx211,AA,AB,AC,JA,JB,JC,1e-5,uk,pk,S);%refine triplet
    uk=uk+alpha*pk; %update uk
    pause(0.01);
    if abs(alpha)<1e-20
        break;
    end
end
% s.mc=1;
% for n=1:s.N+1%Computeukcorrespondingtodifferents.mctogivesamexk
%     E=ComputeE(xk(1:6,n),s);N=ComputeN(xk(1:6,n),0,s);uk(n,1)=s.B'*(E*xk(4:9,n)-N);
% end
end%functionExample211
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=RHS(x,u,s)
E=ComputeE(x,s);
N=ComputeN(x,u,s);
R=E\N;
end%functionRHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=RHSa(r,x,s)
E=ComputeE(x,s);
A=ComputeA(x,s);
R=-E'\(A'*r+s.Q*x(1:4));
end%functionRHSa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E=ComputeE(x,S)
E = [eye(2,2), zeros(2,2); zeros(2,2), [(S.mc+S.m+2*S.I/S.r^2), -S.m*S.l*cos(x(2));-S.m*S.l*cos(x(2)),4/3*S.m*S.l^2]];
end%functionComputeE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N=ComputeN(x,u,S)
N = [x(3);x(4);2*(S.st*S.r*u-S.k*x(3))/S.r^2 - S.m*S.l*sin(x(2))*x(4)^2;S.m*S.g*S.l*sin(x(2))];
end%functionComputeN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=ComputeA(x,S)
A = [zeros(2,2), eye(2,2);[0,-S.m*S.l*(x(6)*sin(x(2))+x(4)^2*cos(x(2)));0,S.m*S.l*(S.g*cos(x(2))-x(5)*sin(x(2)))],[-2*S.k/S.r^2,-2*S.m*S.l*x(4)*sin(x(2));0,0]];
end%functionComputeA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J=ComputeJEx211(utrial,s)
x=s.x0;
u=utrial(1);
J=0.25*s.h*(x'*s.Q*x+u'*s.R*u);
c=.5;
for n=1:s.N
    u=utrial(n);
    if n==s.N
        c=.25;
    end
    f1=RHS(x,u,s);
    f2=RHS(x+s.h*f1/2,u,s);
    f3=RHS(x+s.h*f2/2,u,s);
    f4=RHS(x+s.h*f3,u,s);
    x=x+s.h*(f1/6+(f2+f3)/3+f4/6);
    J=J+c*s.h*(x'*s.Q*x+u'*s.R*u);
end
E=ComputeE(x,s);
J=J+0.5*(x'*E'*s.QT*E*x);
end%functionComputeJEx211

function [AA,AB,AC,JA,JB,JC]= Bracket(ComputeJ,AA,AB,JA,X,P,V)%NumericalRenaissanceCodebase1.0
JB=ComputeJ(X+AB*P,V);
if JB>JA
    %[AA,AB]=Swap(AA,AB);
    temp = AA;
    AA = AB;
    AB = temp;
    %[JA,JB]=Swap(JA,JB);
    temp = JA;
    JA = JB;
    JB = temp;
end
AC=AB+2*(AB-AA);
JC=ComputeJ(X+AC*P,V);
end%functionBracket

function[AB,JB]=Brent(ComputeJ,AA,AB,AC,JA,JB,JC,TOL,X,P,V)
%NumericalRenaissanceCodebase1.0
%INPUT:{AA,AB,AC}bracketaminimumofJ(A)=ComputeJ(X+A*P),withvalues{JA,JB,JC}.
%OUTPUT:ABlocallyminimizesJ(A),withaccuarcyTOL*abs(AB)andvalueJB.
AINC=0;
AL=min(AC,AA);
AR=max(AC,AA);
if(abs(AB-AA)>abs(AC-AB))
    [AA,AC]=Swap(AA,AC);
    [JA,JC]=Swap(JA,JC);
end
for ITER=1:55
    if ITER<3
        AINT=2*(AR-AL);
    end
    TOL1=TOL*abs(AB)+1E-25;
    TOL2=2*TOL1;
    FLAG=0;%Initialize
    AM=(AL+AR)/2;
    if(AR-AL)/2+abs(AB-AM)<TOL2
        ITER
        return;
    end%Checkconvergence
    if(abs(AINT)>TOL1||ITER<3)
    %Performaparabolicfitbasedonpoints{AA,AB,AC}[see(15.2)]
        T=(AB-AA)*(JB-JC);
        D=(AB-AC)*(JB-JA);
        N=(AB-AC)*D-(AB-AA)*T;
        D=2*(T-D);
        if D<0
            N=-N;
            D=-D;
        end
        T=AINT;
        AINT=AINC;
        if(abs(N)<abs(D*T/2)&&N>D*(AL-AB)&&N<D*(AR-AB))
            %AINC=N/D within reasonable range?
            AINC=N/D;AN=AB+AINC;
            FLAG=1;%Success!AINCisnewincrement.
            if(AN-AL<TOL2||AR-AN<TOL2)
                AINC=abs(TOL1)*sign(AM-AB);
            end%FixifANnearends
        end
    end
    %Ifparabolicfitunsuccessful,dogoldensectionstepbasedonbracket{AL,AB,AR}
    if FLAG==0
        if AB>AM
            AINT=AL-AB;
        else
            AINT=AR-AB;
        end
        AINC=0.381966*AINT;
    end
    if abs(AINC)>TOL1
        AN=AB+AINC;
    else
        AN=AB+abs(TOL1)*sign(AINC);
    end
    JN=ComputeJ(X+AN*P,V);
    if JN<=JB%Keep6(notnecessarilydistinct)points
        if AN>AB
            AL=AB;
        else
            AR=AB;
        end%definingtheintervalfromoneiteration
        AC=AA;
        JC=JA;
        AA=AB;
        JA=JB;
        AB=AN;
        JB=JN;%tothenext:
    else%{AL,AB,AR}brackettheminimum
        if AN<AB
            AL=AN;
        else
            AR=AN;
        end%AB=Lowestpoint,mostrecentiftiedw/AA
        if(JN<=JA||AA==AB)%AA=Second-to-lowestpoint.
            AC=AA;
            JC=JA;
            AA=AN;
            JA=JN;%AC=Third-to-lowestpoint
        elseif(JN<=JC||AC==AB||AC==AA)%AN=Newestpoint
            AC=AN;
            JC=JN;%Parabolicfitbasedon{AA,AB,AC}
        end%Goldensectionsearchbasedon{AL,AB,AR}
    end
    %if V
        disp(sprintf('%d%9.5f%9.5f%9.5f%9.5f%9.5f',FLAG,AA,AB,AC,AL,AR));
    %end
end
end%functionBrent
%% 800 RPM Grearteasin Motor Param ID
% 8.66 V battery used as step, output measured in degrees with rotary
% encoder

load(['C:\Users\ethan\OneDrive\Documents\MATLAB\Miscellaneous\InvertedPendulumProj\Arduino\data' '\test7']);

t = data(1,:);
vel = diff(data(3,:))./diff(t);
u = data(2,:);

diffs = find(diff(u) ~= 0); diffs = diffs(1:end-6);
indexes = reshape(diffs,[2,length(diffs)/2]);
B1sum = 0; B2sum = 0; B3sum = 0;

f = @(b,t) -b(1).*exp(b(2).*t)+b(3);
figure; hold on;
plot(t(1:end-1),vel);
plot(t,u);
for i = 1:length(diffs)/2
    ts = t(indexes(1,i):indexes(2,i)-1);
    vels = vel(indexes(1,i):indexes(2,i)-1);
    us = u(indexes(1,i):indexes(2,i)-1);
    
    B{i} = fminsearch(@(b) norm(vels-f(b,ts-ts(1))),[-.5;-20;-.5]);
    
    plot(ts,f(B{i},ts-ts(1)));
    B1sum = B1sum + B{i}(1)/us(2);
    B2sum = B2sum + B{i}(2);
    B3sum = B3sum + B{i}(3)/us(2);
end

B1 = B1sum/(length(diffs)/2);
B2 = B2sum/(length(diffs)/2);
B3 = B3sum/(length(diffs)/2);

G = tf([-B1*B2],[1 -B2]);

te = t(1):.01:t(end);
for i = 1:length(te)
    ue(i) = interp1(t,u,te(i));
end
plot(te,lsim(G,ue,te));


I = pi/2*2700*(7/1000*(25/1000)^4-(6/1000)^4+10/1000*((42/1000)^4-(6/1000)^4));
m = 10/1000+25/1000;
r = 21/1000;
st = -B1*B2*(2*I+m*r^2)/(2*r);
k = B2*(2*I+m*r^2)/(-2);

ui = [-1:.1:-.1, .1:.1:1];
for i = 1:length(ui)
    if ui(i)<0
        vi(i) = min(vel(find(ui(i)-.01 <= u & ui(i)+.01 >= u)));
    else
        vi(i) = max(vel(find(ui(i)-.01 <= u & ui(i)+.01 >= u)));
    end
end
figure;
subplot(2,1,1);
plot(ui,vi);
subplot(2,1,2);
plot(ui,2*(st*r*ui-k*vi)/r^2);
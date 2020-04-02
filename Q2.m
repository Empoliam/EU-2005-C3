clear();

%Constants
global a T d D k N v;

a = 2e-4;
T = 30;
d = 0.3/60;
D = 2;
k = 1/60;
N = 150;
v = 2;

%Maximum time, seconds
tMax = 10000;

initialVal = 0;

%Solve and plot
[X,Y] = ode45(@dL,[0,tMax],initialVal);
plot(X,Y,"k")
xticklabels(round(get(gca,'xtick')./60,0))
xlabel("Time (mins)")
ylabel("Flagellum Length (um)")

%Calculate steady state
qA = (d.* v .* k);
qB = (2 .* D.*d.*k + 2.*D.*k.*N.*v.*a);
qC = ((2.*d.*D.*v) - (2.*D.*k.*N.*v.*a.*T));

steadyState = (-qB + sqrt(qB.^2 - 4.*qA.*qC))./(2.*qA);
disp(steadyState)

%Define ODE
function [out] = dL(~,l)

global a T d D k N v;

J = (k .* N) ./ ( 1 + (k .* l)./v + (k .* l.^2)./(2 .* D) );
out = a .* J .* (T - l) - d;

end

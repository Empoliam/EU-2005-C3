clear();

global a T d D k N v;

a = 2e-4;
T = 30;
d = 0.3/60;
D = 2;
k = 1/60;
N = 150;
v = 2;

tMax = 1000;
dt = 2e-2;
t = 0:dt:tMax;
initialVal = 0.01;

hold on

[X,Y] = solveRK4(@dL,[0,tMax],dt,initialVal);
plot(X,Y)

[W,Z] = ode45(@dL,[0,tMax],initialVal);
plot(W,Z)

hold off

qA = (d.* v);
qB = (2 .* a .* D .* N .* v + 2 .* d .* D);
qC = -2.*a.*D.*N.*v.*T;

steadyState = (-qB + sqrt(qB.^2 - 4.*qA.*qC))/(2.*qA);
disp(steadyState)

function [out] = dL(~,l)

global a T d D k N v;

J = (k .* N) ./ ( ( l .* k ./ v ) + ( (l.^2 .* k) ./ (2 .* D) ) );
out = a .* J .* (T - l) - d;

end

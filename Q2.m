clear();

global a T d D k N v;

a = 2e-4;
T = 30;
d = 0.3/60;
D = 2;
k = 1/60;
N = 150;
v = 2;

tMax = 10000;
dt = 2e-2;

initialVal = 0;

% [t,L] = solveRK4(@dL,[0,tMax],dt,initialVal);
% plot(t,L)

[X,Y] = ode45(@dL,[0,tMax],initialVal);
% v = v/2;
% [W,Z] = ode45(@dL,[tMax,2.*tMax],Y(end));
% plot([X;W],[Y;Z]);
plot(X,Y)

qA = (d.* v .* k);
qB = (2 .* D.*d.*k + 2.*D.*k.*N.*v.*a);
qC = ((2.*d.*D.*v) - (2.*D.*k.*N.*v.*a.*T));

steadyState = (-qB + sqrt(qB.^2 - 4.*qA.*qC))./(2.*qA);
disp(steadyState)

function [out] = dL(~,l)

global a T d D k N v;

%J = (k .* N) ./ ( ( l .* k ./ v ) + ( (l.^2 .* k) ./ (2 .* D) ) );
J = (2.*D.*k.*N.*v)./((2.*D.*(k.*l+v))+ (k.*l.^2.*v));
out = a .* J .* (T - l) - d;

end

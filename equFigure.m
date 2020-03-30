clear();

global a T d D k N v;

a = 2e-4;
T = 30;
d = 0.3/60;
D = 2;
k = 1/60;
N = 150;
v = 2;

l = 0:0.01:15;

plot(l,term1(l),"k")
hold on
plot(l,zeros(1,length(l))+d,"k--")
legend("aJ(T-L)","d")
xlabel("L (um)")
ylabel("Rate (um/s)")
hold off

function [out] = term1(l)

global a T D k N v;

J = (2.*D.*k.*N.*v)./((2.*D.*(k.*l+v))+ (k.*l.^2.*v));
out = a .* J .* (T - l);

end
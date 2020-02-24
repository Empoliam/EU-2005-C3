clear();

b = 50/3;
a = 2e-4;
T = 30;
d = 0.3/60;

dL = @(l) b./l .* a .*(T-l) - d;

tMax = 10000;
dt = 1e-1;
t = 0:dt:tMax;

L = zeros(1,length(t));
L(1) = 1;

i = 1;
while i < length(t)

    L(i+1) = L(i) + dt*dL(L(i));
    i = i + 1;
    
end

plot(t,L)

disp((b.*a.*T)./(d+b.*a))
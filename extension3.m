clear();

%Constants
global a T d D k N v db;

a = 2e-4;
T = 30;
d = 0.3/60;
db = 0.001;
D = 2;
k = 1/60;
N = 150;
v = 2;

%Set up timestep
tMax = 10000;
dt = 2e-2;
t = 0:2e-2:tMax;

%Initialise length arrays
Li = zeros(1,length(t));
Lii = zeros(1,length(t));

%Forward euler solve
iter = 1;
while iter < length(t)/2
   
    Li(iter+1) = Li(iter) + dt.*dLi(Li(iter),Lii(iter));
    Lii(iter+1) = Lii(iter) + dt.*dLii(Li(iter),Lii(iter));
    
    iter = iter+1;
    
end

%Sever flaegellum
Li(iter) = 0;

%Continue solving
while iter < length(t)
   
    Li(iter+1) = Li(iter) + dt.*dLi(Li(iter),Lii(iter));
    Lii(iter+1) = Lii(iter) + dt.*dLii(Li(iter),Lii(iter));
    
    iter = iter+1;
    
end

%Plot
subplot(1,2,1)
plot(t,Li,"k")
xticklabels(round(get(gca,'xtick')./60,0))
xlabel("Time (mins)")
ylabel("Flagellum Length (um)")
title("ODE: Flagellum A")
subplot(1,2,2)
plot(t,Lii,"k")
xticklabels(round(get(gca,'xtick')./60,0))
xlabel("Time (mins)")
ylabel("Flagellum Length (um)")
title("ODE: Flagellum B")

%ODE for flagellum A
function [out] = dLi(li,lii)

global a T d D k N v db;

J = (k .* N) ./ ( 1 + (k .*(li + lii))./v + (k .*(li.^2 + lii.^2))./(2 .* D) );

out = a .* J .* (T - li - lii) - (d+(db * (J * li)./D));

end

%Ode for flagellum B
function [out] = dLii(li,lii)

global a T d D k N v db;

J = (k .* N) ./ ( 1 + (k .*(li + lii))./v + (k .*(li.^2 + lii.^2))./(2 .* D) );

out = a .* J .* (T - li - lii) - (d + (db * (J * lii)./D));

end
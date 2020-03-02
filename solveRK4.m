function [t,X] = solveRK4(F,T,h,init)

tMin = T(1);
tMax = T(2);

t = tMin:h:tMax;
maxIter = length(t);

%Intialize spatial arrays
X = zeros(1,maxIter);
%Intial Values
X(1) = init;

i = 1;
while (i<maxIter)
    %k1
    xKi = F(t(i),X(i));
    %k2
    xKii = F(t(i)+(h./2),X(i)+(xKi./2));
    %k3
    xKj = F(t(i)+(h./2),X(i)+(xKii./2));
    %k4
    xKjj = F(t(i)+h,X(i)+xKj);
    
    X(i+1) = X(i) + (h/6) .* (xKi + 2*xKii + 2*xKj + xKjj);
    
    i = i + 1;
    
end

end
clear()

numX = 1e8;
X = zeros(1,numX);

maxIter = 100;

for i = 1 : maxIter
    
    r = randi([0 1],1,numX);
    r(r==0) = -1;
    
    X = X + r;
    
end

histogram(X)
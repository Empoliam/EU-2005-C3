clear();

a = 2e-4;
T = 30;
d = 0.3/60;
k = 1/60;
N = 150;

vMin = 0.1;
vMax = 10;
vStep = 0.01;
vRange = vMin:vStep:vMax;

DMin = 0.1;
DMax = 80;
DStep =  0.01;
DRange = DMin:DStep:DMax;

[v,D] = meshgrid(vRange,DRange);

qA = (d.* v .* k);
qB = (2 .* D.*d.*k + 2.*D.*k.*N.*v.*a);
qC = ((2.*d.*D.*v) - (2.*D.*k.*N.*v.*a.*T));

steadyState = (-qB + sqrt(qB.^2 - 4.*qA.*qC))./(2.*qA);

imagesc(vRange,DRange,steadyState)
colorbar
set(gca,'YDir','normal')

min(min(steadyState))
max(max(steadyState))
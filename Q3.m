clear();

%Constants
a = 2e-4;
T = 30;
d = 0.3/60;
k = 1/60;
N = 150;

%Range to vary v and D
vMin = 0.1;
vMax = 10;
vStep = 0.01;
vRange = vMin:vStep:vMax;

DMin = 0.1;
DMax = 80;
DStep =  0.01;
DRange = DMin:DStep:DMax;

%Create grid of (v,D)
[v,D] = meshgrid(vRange,DRange);

%Calculate steady state length
qA = (d.* v .* k);
qB = (2 .* D.*d.*k + 2.*D.*k.*N.*v.*a);
qC = ((2.*d.*D.*v) - (2.*D.*k.*N.*v.*a.*T));
steadyState = (-qB + sqrt(qB.^2 - 4.*qA.*qC))./(2.*qA);

%Plot
imagesc(vRange,DRange,steadyState)
colorbar
colormap gray;
title("Phase Diagram for Steady State Flagellum Lengths")
xlabel("v")
ylabel("D")
set(gca,'YDir','normal')
clear()

data = readmatrix("out.csv");

% plot(data(:,1),data(:,2),"k")
% xticklabels(round(get(gca,'xtick')./60,0))
% xlabel("Time (mins)")
% ylabel("Flagellum Length (um)")
 
% subplot(1,2,1)
% plot(data(:,1),data(:,2),"k")
% xticklabels(round(get(gca,'xtick')./60,0))
% xlabel("Time (mins)")
% ylabel("Flagellum Length (um)")
% title("Stochastic: Flagellum A")
% subplot(1,2,2)
% plot(data(:,1),data(:,3),"k")
% xticklabels(round(get(gca,'xtick')./60,0))
% xlabel("Time (mins)")
% ylabel("Flagellum Length (um)")
% title("Stochastic: Flagellum B")

plot(data(:,1),data(:,2),"k")
hold on
plot(data(:,1),data(:,3),"r")
hold off
xlabel("Time (mins)")
ylabel("Flagellum Length (um)")
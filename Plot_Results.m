%% Script to plot results from Johnston, Simpson and Crampin
%% "Predicting population extinction in lattice-based birth-death-movement models"

%% Plot the evolution of the average occupancy from the SSDA PDE, Random Walk and Mean Field ODE
figure(1); hold on; 
plot(0:dt:tEnd,evolutionPDE,'linewidth',2)
plot(0:tEnd,evolutionIBM,'linewidth',2)
plot(odet,odeSolution,'linewidth',2)
legend('SSDA PDE','IBM','ODE')
box on;

%% Plot the evolution of the extinction probability from the SSDA PDE and the Random Walk
figure(2); 
hold on; plot(0:dt:tEnd,extinctionPDE,'linewidth',2)
plot(0:tEnd,extinctionIBM/nRepeats,'linewidth',2)
legend('SSDA PDE','IBM')
box on;

%% Plot the probability density function from the SSDA PDE and the Random Walk
figure(3); hold on;
plot(S,solution,'linewidth',2)
plot(linspace(1/nSites,1,nSites),densityResult(2:end)/nRepeats*nSites,'linewidth',2)
box on;
legend('SSDA PDE','IBM')
xlim([1/nSites 1])
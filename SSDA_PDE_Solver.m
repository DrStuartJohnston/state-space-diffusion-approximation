%% Script to obtain numerical solutions to the SSDA PDE of Johnston, Simpson and Crampin
%% "Predicting population extinction in lattice-based birth-death-movement models"

solution = zeros(nSpaceSites,1);                                            %Initialise PDE solution
spaceVec = linspace(1/nSites,1,numel(solution));                            %Initialise vector of state space
ds = spaceVec(2)-spaceVec(1);                                               %Space between grid points
[~,relIndice] = min(abs(spaceVec-initialDensity));                          %Closest grid point to initial density location
solution(relIndice) = numel(solution);                                      %Set initial condition
IC = solution;                                                              %Initial condition
S = spaceVec;                                                               %Vector of state space

a = zeros(size(solution));                                                  %Initialise subdiagonal coefficients for Thomas_Algorithm
b = zeros(size(solution));                                                  %Initialise diagonal coefficents for Thomas_Algorithm
c = zeros(size(solution));                                                  %Initialise superdiagonal coefficients for Thomas_Algorithm

as = pProliferation_g*nSites^4/((nSites-1)*(nSites-2)*(nSites-3)*(nSites-4))* ...   %Coefficients of a(s)
    (S-1/nSites).*(1-S+1/nSites).*((nSites-2)*(nSites-3)*(nSites-4)/nSites^3 - ...
    (1-S).*(1-S-1/nSites).*(1-S-2/nSites)) + ...
    pProliferation_i*nSites^4/((nSites-1)*(nSites-2)*(nSites-3)*(nSites-4))* ...
    (S-1/nSites).*(1-S+1/nSites).*(1-S).*(1-S-1/nSites).*(1-S-2/nSites);

bs = pDeath_g*nSites^4/((nSites-1)*(nSites-2)*(nSites-3)*(nSites-4))* ...           %Coefficients of b(s)
    (S+1/nSites).*((nSites-1)*(nSites-2)*(nSites-3)*(nSites-4)/nSites^4 - ...
    (1-S-1/nSites).*(1-S-2/nSites).*(1-S-3/nSites).*(1-S-4/nSites)) + ...
    pDeath_i*nSites^4/((nSites-1)*(nSites-2)*(nSites-3)*(nSites-4))* ...
    (S+1/nSites).*(1-S-1/nSites).*(1-S-2/nSites).*(1-S-3/nSites).*(1-S-4/nSites);

source = nSites^4/((nSites-1)*(nSites-2)*(nSites-3)*(nSites-4))*(...                %Coefficients of source term (i.e. if PDE isn't written in conservative form)
    pProliferation_g*((2*S-1-1/nSites)*(1-2/nSites)*(1-3/nSites)*(1-4/nSites) - ...
    (5*S-1-1/nSites).*(1-S).*(1-S-1/nSites).*(1-S-2/nSites)) + ...
    pProliferation_i*(5*S-1-1/nSites).*(1-S).*(1-S-1/nSites).*(1-S-2/nSites) + ...
    pDeath_g*((1-1/nSites)*(1-2/nSites)*(1-3/nSites)*(1-4/nSites) + ...
    (5*S-1+4/nSites).*(1-S-1/nSites).*(1-S-2/nSites).*(1-S-3/nSites)) + ...
    pDeath_i*(1-5*S-4/nSites).*(1-S-1/nSites).*(1-S-2/nSites).*(1-S-3/nSites));

ads = nSites^4/((nSites-1)*(nSites-2)*(nSites-3)*(nSites-4))*(...                   %Coefficients of a'(s)
    pProliferation_g*((1-2*S+2/nSites).*((nSites-2)*(nSites-3)*(nSites-4)/nSites^3 - ...
    (1-S).*(1-S-1/nSites).*(1-S-2/nSites)) + (S-1/nSites).*(1-S+1/nSites).* ...
    ((1-S).*(1-S-1/nSites) + (1-S).*(1-S-2/nSites) + (1-S-1/nSites).*(1-S-2/nSites))) + ...
    pProliferation_i*((1-S+1/nSites).*(1-S).*(1-S-1/nSites).*(1-S-2/nSites) - (S-1/nSites).*(...
    (1-S).*(1-S-1/nSites).*(1-S-2/nSites) + (1-S+1/nSites).*(1-S-1/nSites).*(1-S-2/nSites) + ...
    (1-S+1/nSites).*(1-S).*(1-S-2/nSites) + (1-S+1/nSites).*(1-S).*(1-S-1/nSites))));

bds = nSites^4/((nSites-1)*(nSites-2)*(nSites-3)*(nSites-4))*(...                   %Coefficients of b'(s)
    pDeath_g*((nSites-1)*(nSites-2)*(nSites-3)*(nSites-4)/nSites^4 - ...
    (1-S-1/nSites).*(1-S-2/nSites).*(1-S-3/nSites).*(1-S-4/nSites) + (S+1/nSites).*(...
    (1-S-2/nSites).*(1-S-3/nSites).*(1-S-4/nSites) + (1-S-1/nSites).*(1-S-3/nSites).*(1-S-4/nSites) + ...
    (1-S-1/nSites).*(1-S-2/nSites).*(1-S-4/nSites) + (1-S-1/nSites).*(1-S-2/nSites).*(1-S-3/nSites))) + ...
    pDeath_i*((1-S-1/nSites).*(1-S-2/nSites).*(1-S-3/nSites).*(1-S-4/nSites) - (S+1/nSites).*(...
    (1-S-2/nSites).*(1-S-3/nSites).*(1-S-4/nSites) + (1-S-1/nSites).*(1-S-3/nSites).*(1-S-4/nSites) + ...
    (1-S-1/nSites).*(1-S-2/nSites).*(1-S-4/nSites) + (1-S-1/nSites).*(1-S-2/nSites).*(1-S-3/nSites))));

b(2:end-1) = 1 + 2*D*(as(2:end-1)+bs(2:end-1))*dt/(ds^2) - source(2:end-1)*dt;      %Diagonal coefficients for Thomas Algorithm

b(1) = 1 + 2*D*dt/ds^2*(as(1)+bs(1)) - source(1)*dt + ...                           %Diagonal coefficient for Thomas Algorithm for LHS boundary
    2*ds*(bs(1)-as(1)-D*(ads(1)+bds(1))-pDeath_i/nSites)/(D*(as(1)+bs(1))) * ...
    (-D*dt/ds^2*(as(1)+bs(1)) + (bs(1)-as(1))*dt/(2*ds));

b(end) = 1 + 2*D*dt/ds^2*(as(end)+bs(end)) - source(end)*dt + ...                   %Diagonal coefficient for Thomas Algorithm for RHS boundary
    2*ds*(bs(end)-as(end)-D*(ads(end)+bds(end)))/(D*(as(end)+bs(end))) * ...
    (D*dt/ds^2*(as(end)+bs(end)) + (bs(end)-as(end))*dt/(2*ds));

a(2:end-1) = - D*(as(2:end-1)+bs(2:end-1))*(1)*dt/(ds^2) + dt*(-as(2:end-1)+bs(2:end-1))/(2*ds); %Subdiagonal coefficients for Thomas Algorithm
a(end) = -2*D*dt/ds^2*(as(end)+bs(end));                                                         %Subdiagonal coefficient for Thomas Algorithm for LHS boundary

c(2:end-1) = - D*(as(2:end-1)+bs(2:end-1))*(1)*dt/(ds^2) - dt*(-as(2:end-1)+bs(2:end-1))/(2*ds); %Superdiagonal coefficients for Thomas Algorithm
c(1) = -2*D*dt/ds^2*(as(1)+bs(1));                                                               %Superdiagonal coefficient for Thomas Algorithm for RHS boundary

evolutionPDE = zeros(tEnd/dt+1,1);                                          %Initialise average occupancy obtained from PDE solution
extinctionPDE = evolutionPDE;                                               %Initialise extinction probability
evolutionPDE(1) = mean(spaceVec'.*solution);                                %Initial average occupancy value

count = 1;
    
for t = dt:dt:tEnd                                                          %Iterate over timesteps
    count = count+1;    
    solution = Thomas_Algorithm(a,b,c,solution);                            %Solve Thomas Algorithm to obtain new solution
    evolutionPDE(count) = mean(spaceVec.*solution);                         %New average occupancy
    extinctionPDE(count) = 1 - sum(solution)/numel(solution);               %New extinction probability
end

extinctionPDEFinal = 1 - sum(solution)/numel(solution);                     %Final extinction probability
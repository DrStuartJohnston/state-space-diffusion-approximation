%% Script to obtain numerical solutions to the mean field ODE of Johnston, Simpson and Crampin
%% "Predicting population extinction in lattice-based birth-death-movement models"

ODEParameters.ppg = pProliferation_g;                                       %Parameter structure for grouped proliferation
ODEParameters.ppi = pProliferation_i;                                       %Parameter structure for isolated proliferation
ODEParameters.pdg = pDeath_g;                                               %Parameter structure for grouped death
ODEParameters.pdi = pDeath_i;                                               %Parameter structure for isolated death

[odet,odeSolution] = ode45(@(t,C) def_ODE(t,C,ODEParameters),[0 tEnd],initialDensity_start); %Obtain solution to mean field ODE

function dC = def_ODE(t,C,ODEParameters)                                    %Define mean field ODE

dC = ODEParameters.ppg*C*(1-C)*(1-(1-C)^3) + ...
    ODEParameters.ppi*C*(1-C)^4 - ...
    ODEParameters.pdg*C*(1-(1-C)^4) - ...
    ODEParameters.pdi*C*(1-C)^4;

end
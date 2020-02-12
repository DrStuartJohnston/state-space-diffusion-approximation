%% Highest level script for results presented in Johnston, Simpson and Crampin
%% "Predicting population extinction in lattice-based birth-death-movement models"
%% To perform logistic model simulations ensure pProliferation_i = pProliferation_g
%% and pDeath_i = pDeath_g.

close all
clear all

nSites = 100;                                               %Number of lattice sites
nRepeats = 100;                                            %Number of realisations of the random walk
nSteps = 10000;                                             %Number of time steps in the PDE solution
nSpaceSites = 10000;                                        %Number of grid points in the PDE solution
nDimensions = 2;                                            %Dimension of lattice
pMove = 1;                                                  %Probability of movement
pProliferation_i = 0.03;                                    %Probability of proliferation for isolated individuals
pProliferation_g = 0.01;                                    %Probability of proliferation for grouped individuals
pDeath_i = 0.02;                                            %Probability of death for isolated individuals
pDeath_g = 0.01;                                            %Probability of death for grouped individuals
tEnd = 2000;                                                %Final time in random walk and PDE/ODE solution
dt = tEnd/nSteps;                                           %Time step in PDE solution
D = 1/(2*nSites);                                           %Diffusion coefficient
initialDensity_start = 50/nSites;                           %Proportion of initially occupied lattice sites
initialDensity = initialDensity_start;                      %Proportion of initially occupied lattice sites

if floor(nSites^(1/2)) ~= nSites^(1/2)                      %Ensure lattice is square
    error('Mismatch in dimension: sqrt(nSites) must be an integer')
end

Random_Walk;                                                %Run random walk script
initialDensity = ceil(initialDensity*nSites)/nSites;        %Ensure initial density is a fraction of nSites     
SSDA_PDE_Solver;                                            %Run PDE solver for SSDA PDE
Mean_Field_ODE;                                             %ODE solver for mean field ODE
Plot_Results;                                               %Plot Results

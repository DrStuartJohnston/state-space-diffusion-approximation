%% Script to perform the random walk component of Johnston, Simpson and Crampin
%% "Predicting population extinction in lattice-based birth-death-movement models"

dimensions = nSites^(1/nDimensions)*ones(nDimensions,1);                    %Number of dimensions in random walk
lattice = zeros(nSites,nRepeats);                                           %Initialise lattice

%% Set the initial condition each realisation of the random walk

for j = 1:nRepeats
    lattice(randperm(nSites,ceil(initialDensity*nSites)),j) = 1;            %Set lattice sites where individuals are located
end

evolutionIBM = zeros(tEnd+1,1);                                             %Initialise average occupancy
extinctionIBM = evolutionIBM;                                               %Initialise extinction probability
evolutionIBM(1) = sum(sum(lattice))/(nRepeats*nSites);                      %First average occupancy measurement

%% Perform realisations of the random walk

for j = 1:nRepeats
    if mod(j,round(nRepeats/20))==0
        fprintf('Random walk is %f percent done \n',100*j/nRepeats)     %Progress tracker
    end
    occupiedSites = find(lattice(:,j)==1);                                  %Locations of agents on the lattice
    for i = 1:tEnd                                                          %Iterate over time
        for k = 1:sum(lattice(:,j))                                         %Iterate over population
            pProliferationcheck_i = rand;                                   
            if pProliferationcheck_i < pProliferation_i                     %Check if isolated birth event occurs
                siteNum = ceil(numel(occupiedSites)*rand);                  %Select agent
                selectedSite = occupiedSites(siteNum);                      %Selected agent location (in vector)
                siteCoords = Reverse_Mapping(dimensions,selectedSite);       %Cartesian co-ordinate of agent location
                dircheck = rand;                                            %Proliferation direction
                checkVar = 0;
                %% Check if surrounding sites are occupied (i.e. if agent is isolated or grouped)
                for m = 1:2*nDimensions
                    shift = zeros(nDimensions,1);
                    shift(ceil(m/2)) = (-1)^m;
                    newCoords = siteCoords+shift;
                    newCoords = mod(newCoords-1,dimensions)+1;
                    newSite = Mapping(dimensions,newCoords);                %Maps Cartesian co-ordinates of agent to vector location of agent
                    if lattice(newSite,j) == 1                              %Checks if target site is located
                        checkVar = 1;
                    end
                end
                if checkVar == 0                                            %If agent is isolated
                    %% Check if target site is occupied (i.e. if proliferation event will be unsuccessful)
                    for m = 1:2*nDimensions
                        if dircheck < m/(2*nDimensions) && dircheck > (m-1)/(2*nDimensions)
                            shift = zeros(nDimensions,1);                   %Shift in location for target sites
                            shift(ceil(m/2)) = (-1)^m;
                            newCoords = siteCoords+shift;
                            newCoords = mod(newCoords-1,dimensions)+1;
                            newSite = Mapping(dimensions,newCoords);
                            if lattice(newSite,j) == 0                      %If target site is unoccupied, add new site to agent locations
                                lattice(newSite,j) = 1;
                                occupiedSites = [occupiedSites;newSite];
                                break
                            end
                        end
                    end
                end
            end
            pProliferationcheck_g = rand;
            if pProliferationcheck_g < pProliferation_g                     %Check if grouped proliferation event occurs
                siteNum = ceil(numel(occupiedSites)*rand);                  %Select agent
                selectedSite = occupiedSites(siteNum);                      %Location of agent (in vector)
                siteCoords = Reverse_Mapping(dimensions,selectedSite);      %Cartesian co-ordinate of agent
                dircheck = rand;                                            %Proliferation direction
                checkVar = 2*nDimensions;
                %% Check if surrounding sites are occupied (i.e. if agent is isolated or grouped)
                for m = 1:2*nDimensions
                    shift = zeros(nDimensions,1);
                    shift(ceil(m/2)) = (-1)^m;
                    newCoords = siteCoords+shift;
                    newCoords = mod(newCoords-1,dimensions)+1;
                    newSite = Mapping(dimensions,newCoords);                %Maps Cartesian co-ordinates of agent to vector location of agent
                    if lattice(newSite,j) == 0                              %Checks if target site is located
                        checkVar = checkVar - 1;
                    end
                end
                if checkVar > 0                                             %If agent is grouped
                    %% Check if target site is occupied (i.e. if proliferation event will be unsuccessful)
                    for m = 1:2*nDimensions
                        if dircheck < m/(2*nDimensions) && dircheck > (m-1)/(2*nDimensions)
                            shift = zeros(nDimensions,1);
                            shift(ceil(m/2)) = (-1)^m;
                            newCoords = siteCoords+shift;
                            newCoords = mod(newCoords-1,dimensions)+1;
                            newSite = Mapping(dimensions,newCoords);
                            if lattice(newSite,j) == 0                      %Check if target site is occupied
                                lattice(newSite,j) = 1;
                                occupiedSites = [occupiedSites;newSite];    %If target site is unoccupied, add new site to agent locations 
                                break
                            end
                        end
                    end
                end
            end
            pMovecheck = rand;
            if pMovecheck < pMove                                           %Check if movement event occurs
                siteNum = ceil(numel(occupiedSites)*rand);                  %Select agent
                selectedSite = occupiedSites(siteNum);                      %Location of agent (in vector)
                siteCoords = Reverse_Mapping(dimensions,selectedSite);      %Cartesian co-ordinate of agent
                dircheck = rand;                                            %Movement direction
                %% Check if target site is occupied (i.e. movement event will be unsuccessful)
                for m = 1:2*nDimensions
                    if dircheck < m/(2*nDimensions) && dircheck > (m-1)/(2*nDimensions)
                        shift = zeros(nDimensions,1);
                        shift(ceil(m/2)) = (-1)^m;
                        newCoords = siteCoords+shift;
                        newCoords = mod(newCoords-1,dimensions)+1;
                        newSite = Mapping(dimensions,newCoords);                
                        if lattice(newSite,j) == 0                          %Check if target site is unoccupied
                            lattice(newSite,j) = 1;                         %Change target site to occupied
                            lattice(selectedSite,j) = 0;                    %Change previous site to unoccupied
                            occupiedSites(siteNum) = newSite;               %Update location of agent
                            break
                        end
                    end
                end
            end
            pDeathcheck_i = rand;
            if pDeathcheck_i < pDeath_i                                     %Check if isolated death evemt occurs
                siteNum = ceil(numel(occupiedSites)*rand);                  %Select agent
                selectedSite = occupiedSites(siteNum);                      %Location of agent (in vector) 
                siteCoords = Reverse_Mapping(dimensions,selectedSite);      %Cartesian co-ordinate of agent
                checkVar = 0;
                %% Check if agent is isolated or grouped
                for m = 1:2*nDimensions
                    shift = zeros(nDimensions,1);
                    shift(ceil(m/2)) = (-1)^m;
                    newCoords = siteCoords+shift;
                    newCoords = mod(newCoords-1,dimensions)+1;
                    newSite = Mapping(dimensions,newCoords);
                    if lattice(newSite,j) == 1                              %If one neighbouring site is occupied then agent is not isolated
                        checkVar = 1;
                    end
                end
                if checkVar == 0                                            %If agent is isolated
                    siteNum = ceil(numel(occupiedSites)*rand);
                    selectedSite = occupiedSites(siteNum);
                    lattice(selectedSite,j) = 0;                            %Set lattice site to unoccupied
                    occupiedSites(siteNum) = [];                            %Remove agent location
                end
            end
            if numel(occupiedSites) == 0
                break
            end
            pDeathcheck_g = rand;
            if pDeathcheck_g < pDeath_g                                     %Check if grouped death event occurs
                siteNum = ceil(numel(occupiedSites)*rand);                  %Select agent
                selectedSite = occupiedSites(siteNum);                      %Location of agent (in vector)
                siteCoords = Reverse_Mapping(dimensions,selectedSite);      %Cartesian co-ordinate of agent
                checkVar = 2*nDimensions;
                %% Check if agent is isolated or grouped
                for m = 1:2*nDimensions
                    shift = zeros(nDimensions,1);
                    shift(ceil(m/2)) = (-1)^m;
                    newCoords = siteCoords+shift;
                    newCoords = mod(newCoords-1,dimensions)+1;
                    newSite = Mapping(dimensions,newCoords);
                    if lattice(newSite,j) == 0
                        checkVar = checkVar - 1;                            %If all neighbouring sites are unoccupied then agent is not grouped
                    end
                end
                if checkVar > 0                                             %If agent is grouped
                    siteNum = ceil(numel(occupiedSites)*rand);
                    selectedSite = occupiedSites(siteNum);
                    lattice(selectedSite,j) = 0;                            %Set lattice site to unoccupied
                    occupiedSites(siteNum) = [];                            %Remove agent location
                end
            end
        end
        evolutionIBM(i+1) = evolutionIBM(i+1) + sum(lattice(:,j));          %Calculate average occupancy at time step
        extinctionIBM(i+1) = extinctionIBM(i+1) + sum(sum(lattice(:,j))==0);%Calculate number of populations that have undergone extinction
    end
end

sumLattice = sum(lattice,1);
densityResult = zeros(nSites+1,1);

for i = 1:nSites+1
    densityResult(i) = sum(sumLattice==(i-1));                              %Calculate probability density function
end

evolutionIBM(2:end) = evolutionIBM(2:end)/(nSites*nRepeats);                %Scale occupancy by number of sites and realisations

extinctionIBMFinal = densityResult(1)/sum(densityResult);                   %Final extinction chance
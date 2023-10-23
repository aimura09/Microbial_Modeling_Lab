% Load the COBRA toolbox and the model for cyanobacteria strain 6803

% model = readCbModel('iJN678.mat');
load('iJN678.mat')
model = iJN678;

% Set up DFBA parameters for simulation
simulationTime = 10;  % Simulation time in hours
timeStep = 0.1;  % Time step for simulation in hours
numTimeSteps = simulationTime / timeStep;

% Define the objective function for maximizing pyr production
objectiveRxn = 'EX_pyr_e';
model = changeObjective(model, objectiveRxn);

% Initialize variables to store results
dlactateProduction = zeros(numTimeSteps, 1);
time = zeros(numTimeSteps, 1);

% Define the substrate reaction (e.g., CO2 and ACtex uptake)
substrateRxns = {'EX_ac_e', 'EX_co2_e'};
initConcentrations = [10, 10]
initBiomass = [0.1];
% Define the exchange reaction names
% excRxnNames = model.rxns;

plotRxns = {'EX_ac_e', 'EX_co2_e', 'EX_pyr_e'}

substrateRxns = cellstr(substrateRxns);
[concentrationMatrix, excRxnNames, timeVec, biomassVec]  = dynamicFBA(model, substrateRxns, initConcentrations, initBiomass, timeStep, numTimeSteps, plotRxns);

disp('concentrationMatrix:')
disp(concentrationMatrix)


disp('excRxnNames:')
disp(excRxnNames)

disp('timeVec:')
disp(timeVec)

disp('biomassVec:')
disp(biomassVec)


% % Perform DFBA simulation to maximize pyr production
% for i = 1:numTimeSteps
%     % Set the time constraints for the simulation
%     timeStart = (i - 1) * timeStep;
%     timeEnd = i * timeStep;
% 
%     % Define the time intervals for the simulation
%     timeIntervals = [timeStart, timeEnd];
% 
%     % Convert substrateRxns and initConcentrations to cell arrays of character vectors
%     substrateRxns = cellstr(substrateRxns);
%     % initConcentrations = cellstr(initConcentrations);
%     % initBiomass = cellstr(initBiomass);
% 
%     % Perform DFBA simulation for the current time interval
% 
% 
%     % Store the D-lactate production and time
%     time(i) = i * timeStep;
% end

% Plot pyr production over time
% figure;
% plot(time, dlactateProduction, 'b-', 'LineWidth', 2);
% xlabel('Time (hours)');
% ylabel('pyr Production (mmol/gDW/h)');
% title('pyr Production over Time');
% grid on;

% genetic algorithm to find the optimial substrate uptake (CO2tex and
% ACtex) that produces the most lactate under the given constraints

clc
clear

%Load Model (Cyano model)
load('iJN678.mat')
model = iJN678;

% using mosek solver
changeCobraSolver('mosek', 'LP', 1, 1);

% genetic algorithm parameters ( to be changed )
population_size = 50;
num_generations = 50;
crossover_rate = 0.8;
mutation_rate = 0.02;

% initial set of [popSize x 2] parameters
% initialized to random num 0-1 ( to be changed )
% these define the bounds for CO2tex and ACtex, respectively
current_population = rand(population_size, 2);

% run genetic algorithm loop
for generation = 1:num_generations
    % evaluate fitness of each member of the population

    % select individuals for crossover and reproduction
    % parents = ...

    % apply crossover and mutation to create new generation
    % children = ...

    % current_population = children;
end

% select parameters from the final population set that produce most lactate
%




% objective function / evaluate fitness function
function max_lactate = maximize_lactate(model, c02_bound, acetate_bound)
    model = changeRxnBounds(model, {'CO2tex'}, c02_bound, 'b');
    model = changeRxnBounds(model, {'ACtex'}, acetate_bound, 'b');
    
    model = changeObjective(model, {'BIOMASS_Ec_SynMixo', 'LDH_D'});
    FBAsolution = optimizeCbModel(model, 'max');
        
    max_lactate = FBAsolution.f;
end
%Aiko: changed EX_photon_e bound from -285 to -30

%This code is made to see how biomass growth and D-lactate production
%changes with CO2 and acetate
clc
clear

%Load Model (Cyano model)
load('iJN678.mat')
model = iJN678;

% Adding metabolites for reactions
model = addMetabolite(model, 'ethe_c',...
                'metName', 'Ethylene',...
                'metFormula', 'C2H4',...
                'Charge', 0, ...
                'csense', 'E');
            
model = addMetabolite(model, 'M02035_c',...
                'metName', 'Guanidine',...
                'metFormula', 'CH5N3',...
                'Charge', 0, ...
                'csense', 'E');
            
model = addMetabolite(model, 'pyr5c_c',...
                'metName', '1-pyrroline',...
                'metFormula', 'C5H7NO2',...
                'Charge', 0, ...
                'csense', 'E');

% Adding reactions
model = addReaction(model, 'EthyleneProd', 'reactionFormula',...
    '3 akg_c + 3 o2_c + arg__L_c -> 2 ethe_c + succ_c + M02035_c + 1 pyr5c_c + 7 co2_c + 3 h2o_c',... 
    'reversible', false); 

model = addReaction(model, 'EX_Ethylene_D[e]', 'metaboliteList', {'ethe_c'},...
    'stoichCoeffList', [-1]);

model = addReaction(model, 'DlactateProd', 'reactionFormula',...
    'pyr_c + nadh_c -> nad_c + lac__D_c',... 
    'reversible', false);

model = addReaction(model, 'EX_Dlactate_D[e]', 'metaboliteList', {'lac__D_c'},...
    'stoichCoeffList', [-1]);

model = addReaction(model, 'EX_Guanidine_D[e]', 'metaboliteList', {'M02035_c'},...
  'stoichCoeffList', [1]);

model = addReaction(model, 'EX_pyrroline_D[e]', 'metaboliteList', {'pyr5c_c'},...
    'stoichCoeffList', [1]);

ethylene_idx = find(strcmp(model.rxnNames, 'EX_Ethylene_D[e]'));
lactate_idx = find(strcmp(model.rxnNames, 'EX_Dlactate_D[e]'));

model = changeObjective(model,{'EX_Ethylene_D[e]', 'EX_Dlactate_D[e]'},[1, 1]);

% Nandini's 7002 constraints code
model = changeRxnBounds(model,'O2tex' ,-15,'l');
model = changeRxnBounds(model,'EX_photon_e' ,-30,'b'); %-30 or -100 acceptable
% model = changeRxnBounds(model,'ATPM',1.3,'b'); % removed as no ATPM found in 6803

% Diffusion reactions
model = changeRxnBounds(model, {'GLCtex'}, 0, 'b'); 
model = changeRxnBounds(model, {'FRUtex'}, 0, 'b');
model = changeRxnBounds(model, {'GLUtex'}, 0, 'b'); 
model = changeRxnBounds(model, {'HCO3tex'}, 0, 'b'); 
model = changeRxnBounds(model, {'ALAtex'}, 0, 'b');
model = changeRxnBounds(model, {'ARGtex'}, 0, 'b');
model = changeRxnBounds(model, {'GLYtex'}, 0, 'b');
model = changeRxnBounds(model, {'GLNtex'}, 0, 'b');
model = changeRxnBounds(model, {'HIStex'}, 0, 'b');
model = changeRxnBounds(model, {'LEUtex'}, 0, 'b');
model = changeRxnBounds(model, {'LYStex'}, 0, 'b');
model = changeRxnBounds(model, {'PROtex'}, 0, 'b');
model = changeRxnBounds(model, {'SERtex'}, 0, 'b');
model = changeRxnBounds(model, {'UREAtex'}, 0, 'b');
model = changeRxnBounds(model, {'GLCGLYCtex'}, 0, 'b');
model = changeRxnBounds(model, {'SUCRtex'}, 0, 'b');
model = changeRxnBounds(model, {'PTRCtex'}, 0, 'b');
model = changeRxnBounds(model, {'SPMDtex'}, 0, 'b');
model = changeRxnBounds(model, {'CYNTtex'}, 0, 'b');
model = changeRxnBounds(model, {'CITtex'}, 0, 'b');
model = changeRxnBounds(model, {'MALtpp'}, 0, 'b');
model = changeRxnBounds(model, {'MALtex'}, 0, 'b');
model = changeRxnBounds(model, {'AKGtex'}, 0, 'b');
model = changeRxnBounds(model, {'SUCCtpp'}, 0, 'b');
model = changeRxnBounds(model, {'SUCCtex'}, 0, 'b');
model = changeRxnBounds(model, {'PYRtex'}, 0, 'b');
model = changeRxnBounds(model, {'FUMtpp'}, 0, 'b');
model = changeRxnBounds(model, {'FUMtex'}, 0, 'b');
model = changeRxnBounds(model, {'LDH_D'}, 0, 'b');
model = changeRxnBounds(model,'EX_ac_e',0,'b'); % acetate exchange reaction constrained to 0 for 2-D
model = changeRxnBounds(model, 'EX_Ethylene_D[e]', 0, 'l'); % Nandini constraint from 7002 code
model = changeRxnBounds(model, 'EX_Dlactate_D[e]', 0, 'l'); % Nandini constraint from 7002 code
model = changeRxnBounds(model,'EX_co2_e',-0.1,'b'); % add constraint for 1-dimensional carbon input, starting at first const

num_var = length(model.c);
rng default

options = optimoptions('ga','ConstraintTolerance',1e-4, 'UseParallel',true, 'UseVectorized', true,...
                       'CrossoverFraction',0.85,'PopulationSize',300, 'StallGen', 12,'Generations',5000);

c = model.c; % Objective Coeficients
Aeq = (model.S);
beq = (model.b)';
func = @(v) -objective_infile(v, c);
tic;
% Call Genetic Algorithm
[FBAsolution, fval]  = ga(func,num_var,[],[],Aeq,beq,[],[],[], options);
toc;

objective= -fval;
ethylene=FBAsolution(ethylene_idx);
dlactate=FBAsolution(lactate_idx);

disp(ethylene)
disp(dlactate)
disp(FBAsolution);

function f = objective_infile(v, c)
    f = zeros(1,size(v,1));
    for i = 1:size(v,1)
        f(i) = v(i,:)*c;
    end

end

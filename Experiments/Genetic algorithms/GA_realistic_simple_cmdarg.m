%This code is made to see how biomass growth and D-lactate production
%changes with CO2 and acetate
function GA_realistic_simple_cmdarg(args1,args2,args3,args4,args5,args6)
output_folder = args1; 
ConstraintTolerance = str2num(args2); 
CrossoverFraction = str2num(args3);
PopulationSize = str2num(args4);
StallGen = str2num(args5);
Generations = str2num(args6);

%Load Model (Cyano model)
load('iJN678.mat')
model = iJN678;

options = optimoptions('ga','ConstraintTolerance',ConstraintTolerance, 'UseParallel',true, 'UseVectorized', true,...
                       'CrossoverFraction',CrossoverFraction,'PopulationSize',PopulationSize, ...
                       'StallGen', StallGen,'Generations',Generations);

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

model = changeObjective(model,{'EX_Ethylene_D[e]','EX_Dlactate_D[e]'},[1,1]);

% Nandini's 7002 constraints code
model = changeRxnBounds(model,'O2t' ,-15,'l');
model = changeRxnBounds(model,'EX_photon_e' ,-30,'b'); %-30, -100
model = changeRxnBounds(model,'ATPM',1.3,'b');

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

num_var = length(model.c);
rng default

x = 0; % Defining 'x' value that will be used to plot results (x-axis)
for i = 0.1:0.25:1 % First 'for' loop will model 0 to 20 mmol gDw^{-1} hr^{-1} of flux. This can be changed to model different fluxes  
    x = x+1; % Advancing 'x' for the first iteration
    y = 0; % Defining 'y' value that will be used to plot results (y-axis)
    for j = 0.1:0.25:1 % Second 'for' loop will also model 0 to 20 mmol gDw^{-1} hr^{-1} of flux  
        disp([x, y]);
        y = y+1; % Advancing 'y' for the first iteration

        model = changeRxnBounds(model,'EX_co2_e',-i,'b'); % CO2 exchange reaction will be modeled using the first 'for' loop (i values)
        model = changeRxnBounds(model,'EX_ac_e',-j,'b'); % acetate exchange reaction will be modeled using the second 'for' loop (j values)

        c = model.c; % Objective Coeficients
        Aeq = (model.S);
        beq = (model.b)';
        func = @(v) -objective_infile(v, c);
        tic;
        % Call Genetic Algorithm
        [FBAsolution, fval]  = ga(func,num_var,[],[],Aeq,beq,[],[],[], options);
        toc;
        biomass(x,y)=FBAsolution(687);
        lactate(x,y)=FBAsolution(323);
        %objective(x,y)=-fval;
        objective(x,y)=fval;

        CO2(x)= i; % CO2 is a function of 'x' and will model the i values
        Acetate(y)= j; % acetate is a function of 'y' and will model the j values
        %ethylene(y,x)=FBAsolution(865);
        ethylene(y,x)=FBAsolution(864);
        dlactate(y,x)=FBAsolution(866);
    end
end

subplot(2,2,1)
surf(CO2,Acetate,biomass) % 3-D modeling biomass growth per CO2 and acetate injection
title ('Mixotrophic Biomass Production')
xlabel('CO2 uptake rate)','fontweight','bold','fontsize',11)
ylabel('Acetate uptake rate)','fontweight','bold','fontsize',11)
zlabel('Biomass production)','fontweight','bold','fontsize',11)
axis tight;

subplot(2,2,2)
surf(CO2,Acetate,ethylene) % 3-D modeling biomass growth per CO2 and acetate injection
title ('Ethylene Production')
xlabel('CO2 uptake rate','fontweight','bold','fontsize',11)
ylabel('Acetate uptake rate','fontweight','bold','fontsize',11)
zlabel('Ethylene production','fontweight','bold','fontsize',11)
axis tight;

subplot(2,2,3)
surf(CO2,Acetate,dlactate) % 3-D modeling biomass growth per CO2 and acetate injection
title ('D-lactate Production')
xlabel('CO2 uptake rate','fontweight','bold','fontsize',11)
ylabel('Acetate uptake rate','fontweight','bold','fontsize',11)
zlabel('D-lactate production','fontweight','bold','fontsize',11)
axis tight;

subplot(2,2,4)
surf(CO2,Acetate,objective) % 3-D modelling of the obejctive function value$
title ('objective function value')
xlabel('CO2 uptake rate','fontweight','bold','fontsize', 11)
ylabel('Acetate uptake rate','fontweight','bold','fontsize',11)
zlabel('Biomass production','fontweight','bold','fontsize',11)
axis tight;

folder_path = ['output/' output_folder];
filename = [folder_path '/output_plot_copy.fig'];
saveas(gcf, filename);
filename = [folder_path '/output_plot_copy.png'];
saveas(gcf, filename);

filename = [folder_path '/CO2.mat'];
save(filename, 'CO2')
filename = [folder_path '/Acetate.mat'];
save(filename, 'Acetate')
filename = [folder_path '/biomass.mat'];
save(filename, 'biomass')
filename = [folder_path '/lactate.mat'];
save(filename, 'lactate')
filename = [folder_path '/objective.mat'];
save(filename, 'objective')
filename = [folder_path '/dlactate.mat'];
save(filename, 'dlactate')
filename = [folder_path '/ethylene.mat'];
save(filename, "ethylene")
end

function f = objective_infile(v, c)
    f = zeros(1,size(v,1));
    for i = 1:size(v,1)
        f(i) = v(i,:)*c;
    end

end

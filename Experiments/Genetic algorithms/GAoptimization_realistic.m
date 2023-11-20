%This code is made to see how biomass growth and D-lactate production
%changes with CO2 and acetate
clc
clear all

%Load Model (Cyano model)
load('iJN678.mat')
model = iJN678;
% model = changeObjective(model, {'BIOMASS_Ec_SynMixo','LDH_D'});

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

%constrain biomass
%o2 and photon flux bounded 'b'
%compare fbasolution.f values
%biomass is at 0.08
%d-lactate is at 1
%ethylene is at 0

model = changeObjective(model,{'BIOMASS_Ec_SynMixo','EX_Ethylene_D[e]','EX_Dlactate_D[e]'},[50,10,1]); % Biomass max 53 for d-lactate production for 53,1,10, happens at 13.6
model = changeRxnBounds(model, {'EX_photon_e'},-30,'b'); % Lowest bound is 23 before error
%model = changeRxnBounds(model, {'HCO3tex'},-100,'b'); % Lowest bound is 23 before error


% Pyruvate knockout reactions
model = changeRxnBounds(model, {'PDHcr'}, 0, 'b'); %'Pyruvate dehydrogenase (dihydrolipoamide dehydrogenase) reversible'
model = changeRxnBounds(model, {'PDHbr'}, 0, 'b'); %'Pyruvate dehydrogenase (dihydrolipoamide) reversible'
model = changeRxnBounds(model, {'PDHa'}, 0, 'b'); %'Pyruvate dehydrogenase (lipoamide)'
model = changeRxnBounds(model, {'PDH'}, 0, 'b'); %'Pyruvate dehydrogenase'
model = changeRxnBounds(model, {'PYK3'}, 0, 'b'); %'Pyruvate kinase 3'
model = changeRxnBounds(model, {'PYK5'}, 0, 'b'); %'Pyruvate kinase 5' 
% Pyruvate knockout reactions

% Diffusion reactions
%model = changeRxnBounds(model, {'HCO3E'}, 0, 'b'); 
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
% num_var = 863;
rng default

options = optimoptions('ga','ConstraintTolerance',1e-2,...
                       'CrossoverFraction',0.7,'PopulationSize',3, 'StallGen', 12,'Generations',3);
                        % 'PlotFcn', @gaplotbestf, 'CrossoverFrac',0.2,'PopulationSize',5,'Generations',1);
                        % 'CrossoverFrac',0.5,'PopulationSize',100,'StallGen',125,'Generations',100);
%options = optimoptions('ga','PlotFcns', @gaplotbestfun, 'PlotInterval', 10, 'PopInitRange', [-10 ; 10]);

x = 0; % Defining 'x' value that will be used to plot results (x-axis)
for i = 0.1:0.1:1 % First 'for' loop will model 0 to 20 mmol gDw^{-1} hr^{-1} of flux. This can be changed to model different fluxes  
    x = x+1; % Advancing 'x' for the first iteration
    y = 0; % Defining 'y' value that will be used to plot results (y-axis)
    for j = 0.1:0.1:1 % Second 'for' loop will also model 0 to 20 mmol gDw^{-1} hr^{-1} of flux     
        y = y+1; % Advancing 'y' for the first iteration

        model = changeRxnBounds(model,'EX_co2_e',-i,'b'); % CO2 exchange reaction will be modeled using the first 'for' loop (i values)
        model = changeRxnBounds(model,'EX_ac_e',-j,'b'); % acetate exchange reaction will be modeled using the second 'for' loop (j values)


        c = model.c; % Objective Coeficients
        Aeq = (model.S);
        beq = (model.b)';
        func = @(v) -objective_infile(v, c);

        % Call Genetic Algorithm
        FBAsolution = ga(func,num_var,[],[],Aeq,beq,[],[],[], options);

        biomass(x,y)=FBAsolution(687);
        lactate(x,y)=FBAsolution(323);
        objective(x,y)=FBAsolution(1);

        CO2(x)= i; % CO2 is a function of 'x' and will model the i values
        Acetate(y)= j; % acetate is a function of 'y' and will model the j values
        ethylene(y,x)=FBAsolution(865);
        dlactate(y,x)=FBAsolution(866);
    end
end

subplot(1,3,1)
surf(CO2,Acetate,biomass) % 3-D modeling biomass growth per CO2 and acetate injection
title ('Mixotrophic Biomass Production')
xlabel('CO2 uptake rate (mmol/gDW^{-1}/hr^{-1})','fontweight','bold','fontsize',11)
ylabel('Acetate uptake rate (mmol/gDW^{-1}/hr^{-1})','fontweight','bold','fontsize',11)
zlabel('Biomass production (mmol/gDW^{-1}/hr^{-1})','fontweight','bold','fontsize',11)
axis tight;

subplot(1,3,2)
surf(CO2,Acetate,ethylene) % 3-D modeling biomass growth per CO2 and acetate injection
title ('Ethylene Production')
xlabel('CO2 uptake rate (mmol/gDW^{-1}/hr^{-1})','fontweight','bold','fontsize',11)
ylabel('Acetate uptake rate (mmol/gDW^{-1}/hr^{-1})','fontweight','bold','fontsize',11)
zlabel('Ethylene production (mmol/gDW^{-1}/hr^{-1})','fontweight','bold','fontsize',11)
axis tight;

subplot(1,3,3)
surf(CO2,Acetate,dlactate) % 3-D modeling biomass growth per CO2 and acetate injection
title ('D-lactate Production')
xlabel('CO2 uptake rate (mmol/gDW^{-1}/hr^{-1})','fontweight','bold','fontsize',11)
ylabel('Acetate uptake rate (mmol/gDW^{-1}/hr^{-1})','fontweight','bold','fontsize',11)
zlabel('D-lactate production (mmol/gDW^{-1}/hr^{-1})','fontweight','bold','fontsize',11)
axis tight;

function f = objective_infile(v, c)
    f = zeros(1,size(v,1));
    for i = 1:size(v,1)
        f(i) = v(i,:)*c;
    end

end


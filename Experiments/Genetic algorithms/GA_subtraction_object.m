%This code is made to see how biomass growth and D-lactate production
%changes with CO2 and acetate
clc
clear

options = optimoptions('ga','ConstraintTolerance',1e-4, 'UseParallel',true, 'UseVectorized', true, ... %
                       'CrossoverFraction',0.85,'PopulationSize',500, 'StallGen', 12,'Generations',1000);

load('iJN678.mat')
model = iJN678;

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

model = addReaction(model, 'EthyleneProd', 'reactionFormula',...
    '3 akg_c + 3 o2_c + arg__L_c -> 2 ethe_c + succ_c + M02035_c + 1 pyr5c_c + 7 co2_c + 3 h2o_c',... 
    'reversible', false); 

model = addReaction(model, 'EX_Ethylene_D[e]', 'metaboliteList', {'ethe_c'},...
    'stoichCoeffList', [-1]);

model = addReaction(model, 'EX_Dlactate_D[e]', 'metaboliteList', {'lac__D_c'},...
    'stoichCoeffList', [-1]);

model = addReaction(model, 'EX_Guanidine_D[e]', 'metaboliteList', {'M02035_c'},...
  'stoichCoeffList', [-1]);

model = addReaction(model, 'EX_pyrroline_D[e]', 'metaboliteList', {'pyr5c_c'},...
    'stoichCoeffList', [-1]);


% Nandini's 7002 constraints code
model = changeRxnBounds(model,'O2tex' ,-15,'l');
model = changeRxnBounds(model,'EX_photon_e' ,-30,'b'); %-30, -100

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
model = changeRxnBounds(model, {'EX_Ethylene_D[e]'}, 0, 'l');
model = changeRxnBounds(model, {'EX_Dlactate_D[e]'}, 0, 'l');

% change objective function to the new pseudo lactate reaction
model = changeObjective(model,{'EX_Dlactate_D[e]'});
%model = changeObjective(model,{'LDH_D'})

num_var = length(model.c);
rng default

max_fval = -10000000;
max_ex_ethy = 0;
max_ex_dlac = 0;
max_fval_flux = [];
max_ex_ethy_flux = [];
max_ex_dlac_flux = [];

x = 0;
for i = 0.1:0.25:1
    x = x+1;
    y = 0; 
    for j = 0.1:0.25:1
        y = y+1;
        disp([x, y]);
        disp([i, j])

        model = changeRxnBounds(model,'EX_co2_e',-i,'b');
        model = changeRxnBounds(model,'EX_ac_e',-j,'b');

        c = model.c;
        Aeq = (model.S);
        beq = (model.b)';
        lb = (model.lb);
        ub = (model.ub);
        func = @(v) objective_infile(v, c, i, j);
        tic;
        % Call Genetic Algorithm
        [FBAsolution, fval]  = ga(func,num_var,[],[],Aeq,beq,lb,ub,[], options);
        toc;

        CO2(x)= i;
        Acetate(y)= j;
        
        ethylene(y,x)=FBAsolution(865); % use 'EX_Ethylene_D[e]'
        dlactate(y,x)=FBAsolution(866); % use 'EX_Dlactate_D[e]'
        biomass(y,x)=FBAsolution(687); % use 'Mixotrophic Biomass Ecuation'
        lactate(y,x)=FBAsolution(323); % use 'D-lactate dehydrogenase'
        objective(y,x) = fval; % removed the negative after adding pseudo lactate variable

        if max_fval < fval
            max_fval = fval;
            max_fval_flux = FBAsolution;
            disp('max_fval changed!')
            disp('CO2 Acetate')
            disp([i j])
        end

         if max_ex_ethy < ethylene(y,x)
            max_ex_ethy = ethylene(y,x);
            max_ex_ethy_flux = FBAsolution;
            disp('max_ex_ethy changed!')
            disp('CO2 Acetate')
            disp([i j])
        end

        % lactate increase found
        if max_ex_dlac < dlactate(y,x)
            max_ex_dlac = dlactate(y,x);
            max_ex_dlac_flux = FBAsolution;
            disp('max_ex_dlac changed!')
            disp('CO2 Acetate')
            disp([i j])
        end
        disp('--------------------------')
    end
end

disp("maximum EX_Dlactate_D[e] value detected:")
disp(max_ex_dlac)

output_folder = 'sub_object';
folder_path = ['output/' output_folder];
% filename = [folder_path '/output_plot_copy.fig'];
% saveas(gcf, filename);
% filename = [folder_path '/output_plot_copy.png'];
% saveas(gcf, filename);

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

filename = [folder_path '/model.mat'];
save(filename, 'model')

filename = [folder_path '/max_fval_flux.mat'];
save(filename, 'max_fval_flux')
filename = [folder_path '/max_ex_ethy_flux.mat'];
save(filename, 'max_ex_ethy_flux')
filename = [folder_path '/max_ex_dlac_flux.mat'];
save(filename, 'max_ex_dlac_flux')



% updated to account for pseudo lactate
function f = objective_infile(v, c, co2, ac)
    f = zeros(1,size(v,1));
    for i = 1:size(v,1)
        % calculate the objective function value
        f(i) = co2 + ac - v(i,:)*c;
        %f(i) = co2 + ac - v(i,866);
    end
end

% 
% function f = objective_infile(v, c)
%     disp(size(v))
%     disp(size(c))
%     disp(c)
%     f = zeros(1,size(v,1));
%     for i = 1:size(v,1)
%         f(i) = v(i,:)*c;
%     end
% end
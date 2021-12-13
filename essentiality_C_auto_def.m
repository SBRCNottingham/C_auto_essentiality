%% Gene essentiality analysis for C. autoethanogenum
% Claudio Tomi-Andrino (2021)
% GSM from (10.1049/enb.2018.5003)
initCobraToolbox(false)
backup_model = readCbModel('metaclau.mat');

%% checking GPRs
gpr = backup_model.rxnGeneMat;
gpr_double = full(gpr);

% from the rxn point of view
gpr_list_rxn = num2cell(zeros(size(gpr_double,1),2));
gpr_list_rxn(:,1) = backup_model.rxns;
for i = 1:size(gpr_double,1)
    tmp = num2cell(nonzeros(gpr_double(i,:)));
    tmp_size = size(tmp,1);
    
    if tmp_size > 0
        gpr_list_rxn(i,2) = num2cell(1);
        tmp_indx = find(gpr_double(i,:));
        gpr_list_rxn(i,3) = {num2cell(tmp_indx)};
        tmp_genes = [];
        for j = 1:size(tmp_indx,2)
            tmp_genes = [tmp_genes; backup_model.genes(tmp_indx(1,j))];
        end
        gpr_list_rxn(i,4) = {tmp_genes};
    else
        gpr_list_rxn(i,2) = num2cell(0);
    end 
    clear tmp
    clear tmp_size
    clear tmp_indx
end
clear i
clear tmp_genes
rxn_with_gpr = nnz(cell2mat(gpr_list_rxn(:,2)));
rxn_with_gpr_percent = (rxn_with_gpr/size(backup_model.rxns,1))*100;

% from the gene point of view
gpr_list_gene = num2cell(zeros(size(gpr_double,2),2));
gpr_list_gene(:,1) = backup_model.genes;
for i = 1:size(gpr_double,2)
    tmp = num2cell(nonzeros(gpr_double(:,i)));
    tmp_size = size(tmp,1);
    
    if tmp_size > 0
        gpr_list_gene(i,2) = num2cell(1);
    else
        gpr_list_gene(i,2) = num2cell(0);
    end 
    clear tmp
    clear tmp_size
end
clear i
gene_with_gpr = nnz(cell2mat(gpr_list_gene(:,2)));
gene_with_gpr_percent = (gene_with_gpr/size(backup_model.genes,1))*100;

%% growth on minimal media + CO
% shut down everything except for CO (rxn number 789)
model = backup_model;
model.lb(484) = 0;    % CO2 uptake
model.ub(484) = 0;
model.lb(792) = 0;    % EX_BETA-D-FRUCTOSE
model.ub(792) = 0;

FBA = optimizeCbModel(model);
essential_genes_min_CO = model.genes;

%%% Essential reactions
num_rxns = size(model.rxns,1);
ko_w_feasible_solution_FBA = num2cell(zeros(num_rxns,3));
ko_w_feasible_solution_FBA(:,1) = num2cell(1:num_rxns)';
ko_w_feasible_solution_FBA(:,2) = model.rxns;

for i = 1:num_rxns
    mymodel = model;
    try
        mymodel.lb(i) = 0;
        mymodel.ub(i) = 0;
        solFBA = optimizeCbModel(mymodel);  % solFBA refers to these solutions, whereas FBA is the overall one
        ko_w_feasible_solution_FBA(i,3) = num2cell(solFBA.f);     
    catch
        ko_w_feasible_solution_FBA(i,3) = num2cell(-1);       % infeasible conditions will be flagged with -1
    end
    
    clear mymodel
end
clear i

% select a cut-off --> no flux in the biomass reaction
for i = 1:num_rxns
    if cell2mat(ko_w_feasible_solution_FBA(i,3)) <= 1e-16 
        ko_w_feasible_solution_FBA(i,4) = num2cell(1);      % flagging with '1' 
    else
        ko_w_feasible_solution_FBA(i,4) = num2cell(0); 
    end    
end
clear i

list_essential_rxns_CO = [];
for i =1:size(ko_w_feasible_solution_FBA,1)
    if cell2mat(ko_w_feasible_solution_FBA(i,4)) == 1
        list_essential_rxns_CO = [list_essential_rxns_CO; ko_w_feasible_solution_FBA(i,:)];
    end    
end
clear i

% Add GPR data
for i = 1:size(list_essential_rxns_CO,1)
    tmp_ind = cell2mat(list_essential_rxns_CO(i,1));
    list_essential_rxns_CO(i,5) = gpr_list_rxn(tmp_ind,4);   
end

% pool all genes codifying for an essential enzyme (i.e. appearing in
% column 5). There may be repeated, but it does not matter
list_pooled_essential_genes_CO = [];
for i = 1:size(list_essential_rxns_CO,1)
    if isempty(list_essential_rxns_CO{i,5}) == 0
    list_pooled_essential_genes_CO = [list_pooled_essential_genes_CO; list_essential_rxns_CO{i,5}];
    end
end

for i = 1:size(essential_genes_min_CO,1)
    tmp_i = string(essential_genes_min_CO(i,1));
    for j = 1:size(list_pooled_essential_genes_CO,1)
        tmp_j = string(list_pooled_essential_genes_CO(j,1));
        if strcmp(tmp_i,tmp_j) == 1
            essential_genes_min_CO(i,2) = list_pooled_essential_genes_CO(j,1);
            essential_genes_min_CO(i,3) = num2cell(1);
        end
    end
end

for i = 1:size(essential_genes_min_CO,1)
   if isempty(essential_genes_min_CO{i,3})
       essential_genes_min_CO(i,3) = num2cell(0);
   end    
end

list_genes_min_CO = cell2mat(essential_genes_min_CO(:,3));
clear ko_w_feasible_solution_FBA
clear FBA

%% Gene essentiality comparison (in vivo vs. in silico). Matthew's correlation coefficient
% External files including the in vivo data have been arranged to only
% consider the same genes as the GSM, in the same order.
% some genes have not been tested for all growth conditions, but
% essentiality has been assumed if in the previous growth conditions they
% were essential
essentiality_in_vivo_data = table2cell(readtable('essentiality_in_vivo.txt'));

% MM + CO
trueMat = cell2mat(essentiality_in_vivo_data(:,4));
predictedMat = list_genes_min_CO;
[MCC, TP, TN, FP, FN] = calMCC_mod(trueMat,predictedMat);

MCC_MM_CO = MCC;
clear trueMat
clear predictedMat




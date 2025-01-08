function [geneEssentiality,rxnEssentiality] = essentialityAnalysis(model,resPath,modelName)
% Script to analyse the essential genes and reactions in a model
% Input 
% model  - Model to analyze
% resPath - Path to store the results files
% modelName - Name of the model to store
% Output :
% geneEssentiality - Information about the essential genes and their
% subsystem
% reactionEssentiality - Information about the essential reactions and
% their subsystems
[grR_genes,grKO_genes,grWT_genes] = singleGeneDeletion(model);
[grR_rxns,grKO_rxns,grWT_rxns] = singleRxnDeletion(model);
EssGenes = model.genes(grKO_genes < 0.1* grWT_genes);
EssRxns = model.rxns(grKO_rxns < 0.1* grWT_rxns);

EssRxns(:,2) = model.rxnNames(findRxnIDs(model,EssRxns(:,1)));
EssRxns(:,3) = model.subSystems(findRxnIDs(model,EssRxns(:,1)));
EssRxns = subSystemSeparation(EssRxns);

geneEssentiality = struct();
geneEssentiality.grR = grR_genes;
geneEssentiality.grWT = grWT_genes;
geneEssentiality.grKO = grKO_genes;
geneEssentiality.EssGenes = EssGenes;

rxnEssentiality = struct();
rxnEssentiality.grR = grR_rxns;
rxnEssentiality.grWT = grWT_rxns;
rxnEssentiality.grKO = grKO_rxns;
rxnEssentiality.EssRxns = EssRxns;

writetable(cell2table(geneEssentiality.EssGenes),[resPath,'EssentialityAnalysis.xlsx'],'sheet',['EssG_',modelName]);
writetable(cell2table(rxnEssentiality.EssRxns),[resPath,'EssentialityAnalysis.xlsx'],'sheet',['EssR_',modelName]);

end

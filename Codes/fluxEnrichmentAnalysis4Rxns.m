
% Script to evaluate the flux enrichment analysis for a given reactions 
function feaRxnResultSig = fluxEnrichmentAnalysis4Rxns(enrichRxnSet,model,resPath,sheetName)
% Input
% enrichRxnSet : Reactions to study the flux enrichment analysis
% model : Model to analyse
% resPath : Path to save the results
% sheetName : SheetName to save the results
% Output
% feaRxnResultSig : Signficant enriched subsystems for the provided
% reactions
if length(enrichRxnSet) >= 1
    enrichRxnIDs = findRxnIDs(model,enrichRxnSet(~cellfun('isempty', enrichRxnSet)));
    enrichRxnIDs(find(enrichRxnIDs ==0)) = [];
    feaRxnResult = FEA(model,enrichRxnIDs,'subSystems');
    feaRxnResultSig = feaRxnResult([1; find(cell2mat(feaRxnResult(2:end,2)) < 0.05)+1],3:5);
    feaRxnResultSig(1,4) = {'EnrichmentSize'};
    feaRxnResultSig(2:end,4) = num2cell(cell2mat(feaRxnResultSig(2:end,2))./cell2mat(feaRxnResultSig(2:end,3)));
else
    disp('No reaction is present')
    feaRxnResultSig ={};
end
writetable(cell2table(feaRxnResultSig),[resPath,'FEA_results.xlsx'],'sheet',sheetName);
end
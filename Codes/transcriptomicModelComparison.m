% Script to compare the two GIMME models to explain the active Rxns,
% enrich/repress Rxns (with respect to M1 (M2/M1)), enrich/rrepress RxnSpan and path to save all these file
function modelCompResults = transcriptomicModelComparison(model1,model2,resPath,M1,M2)
%Input 
% model1 : Transcriptomics integrated model 1
% model2 : Transcriptomics integrated model 2
% resPath : Paths to store all the reuslts file including results of FVA,
% enrich, active Rxns
% M1 - Name of the model 1 
% M2 - name of the model 2
% Output
% modelCompResults : Struct with all the analysis results as fields


%Analysis of reaction varied 

%% Active reaction analysis
activeRxn_M1 = setdiff(model1.rxns,model2.rxns);
activeRxn_M2 = setdiff(model2.rxns,model1.rxns);
% Combining rxnNames and Subsystems
activeRxn_M1(:,2) = model1.rxnNames(findRxnIDs(model1,activeRxn_M1(:,1)));
activeRxn_M1(:,3) = model1.subSystems(findRxnIDs(model1,activeRxn_M1(:,1)));
activeRxn_M2(:,2) = model2.rxnNames(findRxnIDs(model2,activeRxn_M2(:,1)));
activeRxn_M2(:,3) = model2.subSystems(findRxnIDs(model2,activeRxn_M2(:,1)));

activeRxn_M1 = subSystemSeparation(activeRxn_M1);
activeRxn_M2 = subSystemSeparation(activeRxn_M2);

%% Growth analysis
sol_M1 = optimizeCbModel(model1);
gr_M1 = sol_M1.f;
[minF_M1,maxF_M1] = fluxVariability(model1);

sol_M2 = optimizeCbModel(model2);
gr_M2 = sol_M2.f;
[minF_M2,maxF_M2] = fluxVariability(model2);

% Analyzing difference in flux variation between two different condition
unionRxns = union(model1.rxns,model2.rxns);
fluxVariation = zeros(length(unionRxns),2);
Fluxspan_M1 = maxF_M1 - minF_M1;
Fluxspan_M2 = maxF_M2 - minF_M2;
fluxSpanVariation = zeros(length(unionRxns),2);

for i = 1:length(unionRxns)
    if findRxnIDs(model1,unionRxns(i)) ~= 0 
    fluxVariation(i,1) = maxF_M1(findRxnIDs(model1,unionRxns(i)));
    fluxSpanVariation(i,1) = Fluxspan_M1(findRxnIDs(model1,unionRxns(i)));
    end
    if findRxnIDs(model2,unionRxns(i)) ~= 0
    fluxVariation(i,2) = maxF_M2(findRxnIDs(model2,unionRxns(i)));
    fluxSpanVariation(i,2) = Fluxspan_M2(findRxnIDs(model2,unionRxns(i)));
    end
end

[enrichRxns,repressRxns] = fluxVariationEstimation(fluxVariation,model1,model2);
[enrichRxnSpan,repressRxnSpan] = fluxVariationEstimation(fluxSpanVariation,model1,model2);

% Exporting the data
writetable(cell2table(activeRxn_M1),[resPath,'FluxCompAnalysis.xlsx'],'sheet',['Active_', M1]);
writetable(cell2table(activeRxn_M2),[resPath,'FluxCompAnalysis.xlsx'],'sheet',['Active_',M2]);
writetable(cell2table(enrichRxns),[resPath,'FluxCompAnalysis.xlsx'],'sheet',['Enriched@',M2]);
writetable(cell2table(repressRxns),[resPath,'FluxCompAnalysis.xlsx'],'sheet',['Repressed@',M2]);
writetable(cell2table(enrichRxnSpan),[resPath,'FluxCompAnalysis.xlsx'],'sheet',['EnrichedSpan@',M2]);
writetable(cell2table(repressRxnSpan),[resPath,'FluxCompAnalysis.xlsx'],'sheet',['RepressedSpan@',M2]);

writetable(cell2table([model1.rxns,num2cell(minF_M1),num2cell(maxF_M1)]),[resPath,'FluxVariabilityAnalysis.xlsx'],'sheet',['FVA_',M1]);
writetable(cell2table([model2.rxns,num2cell(minF_M2),num2cell(maxF_M2)]),[resPath,'FluxVariabilityAnalysis.xlsx'],'sheet',['FVA_',M2]);

modelCompResults = struct();
modelCompResults.activeM1 = activeRxn_M1;
modelCompResults.activeM2 = activeRxn_M2;
modelCompResults.enrichRxns = enrichRxns;
modelCompResults.repressRxns = repressRxns;
modelCompResults.enrichRxnSpan = enrichRxnSpan;
modelCompResults.repressRxnSpan = repressRxnSpan;

end
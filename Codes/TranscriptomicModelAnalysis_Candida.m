% The script explains the reconstruction of transcriptomic data integration
% to the Candida models and simulating the model using media and
% analyze the differently active reactions  

% The analysis can be extended to the mouse-based and HUVEC- cell line
% based studies
%% Creating transcriptomics integrated Candida models
% Loading C.albicans model
CA = readCbModel('../Models/iRV781_annotated.mat');

% Loading the transcriptomic Data
rawData = readtable('../Transcriptomic Data Files/transcriptomic_Candida_OKF6_SC0534.txt');
genes = table2cell(rawData(:,1));

% 5 hours Control - CA - OKF6
Data5hrControl = table2array(rawData(:,6:7));
[thre5hrControl,GIMME5hrControl] = transcriptomicIntegration2model(CA,Data5hrControl,genes,"GIMME");

% 5 hours Infected - CA - OKF6
Data5hrInfected = table2array(rawData(:,8:9));
[thres5hrInfected,GIMME5hrInfected] = transcriptomicIntegration2model(CA,Data5hrInfected,genes,"GIMME");

%% Simulating a integrated model under DMEM media
% DMEM media is used to grow Candida albicans as a negative control in the
% study
media = readtable('../Media/DMEM media.txt');
mediaNames = table2cell(media(:,1));
mediaValue = table2array(media(:,2));

CA_DMEM = changeRxnBounds(CA,mediaNames,mediaValue,'l');

% Creating integrated model with DMEM media
GIMME5hrControl_DMEM = transcriptomicIntegration2model(CA_DMEM,Data5hrControl,genes,"GIMME");
GIMME5hrInfected_DMEM = transcriptomicIntegration2model(CA_DMEM,Data5hrInfected,genes,"GIMME");

% Flux variability analysis
sol_5hControl = optimizeCbModel(GIMME5hrControl);
gr_5hControl = sol_5hControl.f;
[minF_5hControl,maxF_5hControl] = fluxVariability(GIMME5hrControl);

sol_5hInfected = optimizeCbModel(GIMME5hrInfected);
gr_5hInfected = sol_5hInfected.f;
[minF_5hInfected,maxF_5hInfected] = fluxVariability(GIMME5hrInfected);

resPath = '../candida analysis';
model1Name = '5hrControl';
model2Name = '5hrInfected';
GIMME_FVA_OKF6_5hr_comp = transcriptomicModelComparison(GIMME5hrControl_DMEM,GIMME5hrInfected_DMEM,resPath,model1Name,model2Name);

% Essentiality Studies
[essGene_OKF6_5hrControl,essRxn_OKF6_5hrControl] = essentialityAnalysis(GIMME5hrControl_DMEM,resPath,'5hrControl-OKF6');
[essGene_OKF6_5hrInfected,essRxn_OKF6_5hrInfected] = essentialityAnalysis(GIMME5hrInfected_DMEM,resPath,'5hrInfected-OKF6');

% FEA for varying reaction
fea_OKF6_5hr_CvI_C_activeRxns = fluxEnrichmentAnalysis4Rxns(GIMME_FVA_OKF6_5hr_comp.activeM1(:,1),GIMME5hrControl_DMEM,resPath,'OKF6_5hr_CvI_C_active');
fea_OKF6_5hr_CvI_activeRxns = fluxEnrichmentAnalysis4Rxns(GIMME_FVA_OKF6_5hr_comp.activeM2(:,1),GIMME5hrInfected_DMEM,resPath,'OKF6_5hr_CvI_active');
fea_OKF6_5hr_CvI_enrichRxns = fluxEnrichmentAnalysis4Rxns(GIMME_FVA_OKF6_5hr_comp.enrichRxns(:,1),GIMME5hrInfected_DMEM,resPath,'OKF6_5hr_CvI_enrich');
fea_OKF6_5hr_CvI_repressRxns = fluxEnrichmentAnalysis4Rxns(GIMME_FVA_OKF6_5hr_comp.repressRxns(:,1),GIMME5hrInfected_DMEM,resPath,'OKF6_5hr_CvI_repress');
fea_OKF6_5hr_CvI_enrichSpanRxns = fluxEnrichmentAnalysis4Rxns(GIMME_FVA_OKF6_5hr_comp.enrichRxnSpan(:,1),GIMME5hrInfected_DMEM,resPath,'OKF6_5hr_CvI_enrichSpan');
fea_OKF6_5hr_CvI_repressSpanRxns = fluxEnrichmentAnalysis4Rxns(GIMME_FVA_OKF6_5hr_comp.repressRxnSpan(:,1),GIMME5hrInfected_DMEM,resPath,'OKF6_5hr_CvI_repressSpan');

% FEA for essentiality analysis
rescuedRxns_CvI_OKF6_5hr = essRxn_OKF6_5hrControl.EssRxns(find(~ismember(essRxn_OKF6_5hrControl.EssRxns(:,1),essRxn_OKF6_5hrInfected.EssRxns(:,1))),:);
neededRxns_CvI_OKF6_5hr = essRxn_OKF6_5hrInfected.EssRxns(find(~ismember(essRxn_OKF6_5hrInfected.EssRxns(:,1),essRxn_OKF6_5hrControl.EssRxns(:,1))),:);
fea_rescuedRxns_CvI_OKF6_5hr = fluxEnrichmentAnalysis4Rxns(rescuedRxns_CvI_OKF6_5hr(:,1),GIMME5hrInfected_DMEM,resPath,'rescuedRxns_CvI_OKF6_5hr');
fea_neededRxns_CvI_OKF6_5hr = fluxEnrichmentAnalysis4Rxns(neededRxns_CvI_OKF6_5hr(:,1),GIMME5hrInfected_DMEM,resPath,'neededRxns_CvI_OKF6_5hr');

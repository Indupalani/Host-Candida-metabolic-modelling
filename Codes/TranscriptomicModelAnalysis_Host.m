% The script explains the reconstruction of transcriptomic data integration
% to the cell-line based models and simulating the model using DMEM media and
% analyze the differently active reactions  

% The analysis can be extended to the mouse-based and HUVEC- cell line
% based studies
%% Cell line - HUVEC
% load Recon model
load('../Models/Recon3DModel_301.mat');
host = Recon3DModel;

% load the transcriptomics data
rawData = readtable('../Transcriptomic Data Files/transcriptomic_host_HUVEC_SC0534.txt','ReadVariableNames',false);
genes = table2cell(rawData(2:end,1));

% 5 hours Control - host
Data5hrControl = str2double(table2cell(rawData(2:end,2:3)));
[thre5hrControl,GIMME5hrControl] = transcriptomicIntegration2model(host,Data5hrControl,genes,"GIMME");

% 5 hours infected - Host
Data5hrInfected= str2double(table2cell(rawData(2:end,4:5)));
[thres5hrInfected,GIMME5hrInfected] = transcriptomicIntegration2model(host,Data5hrInfected,genes,"GIMME");

% DMEM media is used to grow Candida albicans as a negative control in the
% study
media = readtable('/home/indumathi/Desktop/PhD/Research work/IITM - NIRRH/Transcriptomic data integration/DMEM media.txt');
mediaNames = table2cell(media(:,1));
mediaValue = table2array(media(:,2));
mediaNames = strrep(mediaNames,'(e)','[e]');
host_DMEM = changeRxnBounds(host,host.rxns(findExcRxns(host)),0,'l');
host_DMEM = changeRxnBounds(host_DMEM,mediaNames,mediaValue,'l');

% Model with the applied media
[thre5hrControl_DMEM,GIMME5hrControl_DMEM] = transcriptomicIntegration2model(host_DMEM,Data5hrControl,genes,"GIMME");
[thres5hrInfected_DMEM,GIMME5hrInfected_DMEM] = transcriptomicIntegration2model(host_DMEM,Data5hrInfected,genes,"GIMME");

[thre5hrControl_DMEM_1,GIMME5hrControl_DMEM_1] = transcriptomicIntegration2model(host_DMEM_1,Data5hrControl,genes,"GIMME");
[thres5hrInfected_DMEM_1,GIMME5hrInfected_DMEM_1] = transcriptomicIntegration2model(host_DMEM_1,Data5hrInfected,genes,"GIMME");

% Analysis using transcriptomicModelComparison function
resPath = '../cell line analysis';
model1Name = '5hrControl';
model2Name = '5hrInfected';
GIMME_FVA_HUVEC_5hr_comp = transcriptomicModelComparison(GIMME5hrControl_DMEM,GIMME5hrInfected_DMEM,resPath,model1Name,model2Name);

% Essentiality Studies
[essGene_HUVEC_5hrControl,essRxn_HUVEC_5hrControl] = essentialityAnalysis(GIMME5hrControl_DMEM,resPath,'5hrControl-HUVEC');
[essGene_HUVEC_5hrInfected,essRxn_HUVEC_5hrInfected] = essentialityAnalysis(GIMME5hrInfected_DMEM,resPath,'5hrInfected-HUVEC');

% FEA for varying reaction
fea_HUVEC_5hr_CvI_C_activeRxns = fluxEnrichmentAnalysis4Rxns(GIMME_FVA_HUVEC_5hr_comp.activeM1(:,1),GIMME5hrControl_DMEM,resPath,'HUVEC_5hr_CvI_C_active');
fea_HUVEC_5hr_CvI_activeRxns = fluxEnrichmentAnalysis4Rxns(GIMME_FVA_HUVEC_5hr_comp.activeM2(:,1),GIMME5hrInfected_DMEM,resPath,'HUVEC_5hr_CvI_active');
fea_HUVEC_5hr_CvI_enrichRxns = fluxEnrichmentAnalysis4Rxns(GIMME_FVA_HUVEC_5hr_comp.enrichRxns(:,1),GIMME5hrInfected_DMEM,resPath,'HUVEC_5hr_CvI_enrich');
fea_HUVEC_5hr_CvI_repressRxns = fluxEnrichmentAnalysis4Rxns(GIMME_FVA_HUVEC_5hr_comp.repressRxns(:,1),GIMME5hrInfected_DMEM,resPath,'HUVEC_5hr_CvI_repress');
fea_HUVEC_5hr_CvI_enrichSpanRxns = fluxEnrichmentAnalysis4Rxns(GIMME_FVA_HUVEC_5hr_comp.enrichRxnSpan(:,1),GIMME5hrInfected_DMEM,resPath,'HUVEC_5hr_CvI_enrichSpan');
fea_HUVEC_5hr_CvI_repressSpanRxns = fluxEnrichmentAnalysis4Rxns(GIMME_FVA_HUVEC_5hr_comp.repressRxnSpan(:,1),GIMME5hrInfected_DMEM,resPath,'HUVEC_5hr_CvI_repressSpan');

% FEA for essentiality analysis
rescuedRxns_CvI_HUVEC_5hr = essRxn_HUVEC_5hrControl.EssRxns(find(~ismember(essRxn_HUVEC_5hrControl.EssRxns(:,1),essRxn_HUVEC_5hrInfected.EssRxns(:,1))),:);
neededRxns_CvI_HUVEC_5hr = essRxn_HUVEC_5hrInfected.EssRxns(find(~ismember(essRxn_HUVEC_5hrInfected.EssRxns(:,1),essRxn_HUVEC_5hrControl.EssRxns(:,1))),:);
fea_rescuedRxns_CvI_HUVEC_5hr = fluxEnrichmentAnalysis4Rxns(rescuedRxns_CvI_HUVEC_5hr(:,1),GIMME5hrInfected_DMEM,resPath,'rescuedRxns_CvI_HUVEC_5hr');
fea_neededRxns_CvI_HUVEC_5hr = fluxEnrichmentAnalysis4Rxns(neededRxns_CvI_HUVEC_5hr(:,1),GIMME5hrInfected_DMEM,resPath,'neededRxns_CvI_HUVEC_5hr');


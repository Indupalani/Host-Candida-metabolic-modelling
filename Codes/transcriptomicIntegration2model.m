function [options,integratedModel] = transcriptomicIntegration2model(model,data,genes,method)
% The function integrates the transcriptomic data into the modelusing
% either iMAT or GIMME (0.1*Median as a threshold).

%Mean value of the data
transcriptomicMean = mean(data,2);
transcriptomicMean(isnan(transcriptomicMean)) = 0;

% Assigning to expression values
expression.rawValue = data;
expression.value = transcriptomicMean;
expression.gene = genes;
[expressionRxns parsedGPR] = mapExpressionToReactions(model, expression);

%% iMAT method
if method == "iMAT"
    
% threshold assignment
lbThreshold = 0.05*median(transcriptomicMean);
ubThreshold = 0.15*median(transcriptomicMean);

% Assigning that to solver
options.solver = 'iMAT';
options.threshold_lb = lbThreshold;
options.threshold_ub = ubThreshold;
options.expressionRxns = expressionRxns;

integratedModel = createTissueSpecificModel(model, options);

%% GIMME model
elseif method == "GIMME"
% threshold assignment
threshold = 0.1*median(transcriptomicMean);

% Assigning that to solver
options.solver = 'GIMME';
options.threshold = threshold;
options.expressionRxns = expressionRxns;

integratedModel = createTissueSpecificModel(model, options);
else
    print("Specify the method - iMAT or GIMME");
end
end
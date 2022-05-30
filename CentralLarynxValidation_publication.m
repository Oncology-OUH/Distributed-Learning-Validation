function CentralLarynxValidation
%I have changed the Central script into a function to avoid that it could
%be influenced by variables that might be present in the workspace when it
%is executed. /CaB

% Central algorithm

%A few initializations
addpath('./src');
Init.clients_names = {'Node_server1','Node_server2','Node_server3'}; % list the allowed clients for this analysis
Init.NameCentralMPIClient='catcomCentral1MPIClient';
Init.PortNumberCentralMPIClient='8080';
Init.FileNameSparql_RDFSTore='sparql_queries\larynx_model_160920.sparql';
Init.algorithm = 'stratified_cox';
%Init.FeaturesToOptimizeFrom={'age_start_rt','hemoglobin','eqd2t','T2','T3','T4','Nplus','genderMale','NonGlottic'};
Init.FeaturesToOptimizeFrom={'egelmeer_pi'};
%Init.FeaturesToOptimizeFrom={'age_start_rt','genderMale','T2','T3','T4','Nplus'}; 
Init.FeaturesLevelSets={{'T2','T3','T4'}};
Init.PI_coefficients=table(0.04136,-0.4000,-0.03506,0.1943,0.7965,1.4546,0.3764,0.8353,0.2695,...
  'VariableNames',...
  {'age_start_rt','hemoglobin','eqd2t','T2','T3','T4','Nplus','genderMale','NonGlottic'});
Init.ImputationVariables={'age_start_rt','hemoglobin','eqd2t','t_stage','n_stage','gender','tumour_loc'};
Init.ImputationVariablesTypes={'numeric','numeric','numeric','factor','factor','factor','factor'};
Init.HardCensoringTime=60;
Init.GlobalSeed=21;
Init.NbootCV=2;
Init.NbootModelInfo=2000;
Init.NAllowedMissing=3;
Init.ConvergTollerance=1e-6;
Init.alpha=0.05;
Init.NumerAtRiskTimePoints=0:6:10*12; %Six month interval of number at risk reported values
Init.PI_thresholds = []; %The PI values are used to divide the Kaplan Meier plots of data and model into risk groups.

% Init.PI_thresholds.catcomManData1 = [-1.6782,-0.7554]; %Whole validation
% Init.PI_thresholds.catcomdata1 = [-1.7487,-0.7508]; %Whole validation
% Init.PI_thresholds.catcomManData1 = [-1.6117,-0.7254]; %Curative cohort
% Init.PI_thresholds.catcomdata1 = [-1.8370,-1.0704]; %Curative cohort

Init.TotalDoseRange = [0,100];
Init.TotalFracRange = [0,60];

Init.CalPlotTimePoints = [24,60];
Init.NrCalPlotTimePoints = 6;
Init.FileSaveParamSelect='./model_results/modelCV';
Init.FileSaveInfo='./Model_results/modelInfo';
Init.FigurSavePath='./figures/larynx';

tic;

%The following line of code is a call to the communication server. The aim
%is to inform the communication server that the central code is alive and
%which local nodes should be included in the model optimization. The
%function is waiting until it gets a handshake from the local nodes that
%they are ready for model optimizations. The needed nodes are defined in
%the variable “clients_names”. Thus after a return from this call, all
%needed computers are active and ready for optimization
mpi = ServerConnectionSetup(Init.NameCentralMPIClient,strjoin(Init.clients_names,','),true,Init.PortNumberCentralMPIClient);

model = struct;
model.algorithm = Init.algorithm;
%model.features = Init.FeaturesToOptimizeFrom;
model.PI_coefficients=Init.PI_coefficients;
model.sparql_query = readTextFile(Init.FileNameSparql_RDFSTore);
model.normalise=0;
model.randseed = Init.GlobalSeed;
model.ImputationVariables=Init.ImputationVariables;
model.ImputationVariablesTypes=cell2table(Init.ImputationVariablesTypes,'VariableNames',Init.ImputationVariables);
model.validation='bootstrap';
model.NAllowedMissing=Init.NAllowedMissing;
model.HardCensoringTime=Init.HardCensoringTime;
model.ConvergTollerance=Init.ConvergTollerance;
model.NumerAtRiskTimePoints=Init.NumerAtRiskTimePoints;
model.TotalDoseRange=Init.TotalDoseRange;
model.TotalFracRange=Init.TotalFracRange;

%model.NModel = size(model.varselect,1);
model.beta=[];
varaiables_in = ff2n(length(Init.FeaturesToOptimizeFrom));
varaiables_in = logical(varaiables_in(2:end,:));
counter=0;
for j=1:size(varaiables_in,1)
  varnames=Init.FeaturesToOptimizeFrom(varaiables_in(j,:));
  %Check whether all level pairs are all in or all out
  OkToInclude=false(length(Init.FeaturesLevelSets),1);
  for i=1:length(Init.FeaturesLevelSets)
    temp=setdiff(Init.FeaturesLevelSets{i},varnames);
    if length(temp)==length(Init.FeaturesLevelSets{i}) || isempty(temp)
      %All levels are either all out or in for this specific set of levels
      OkToInclude(i)=true;
    end
  end
  if all(OkToInclude)
    counter=counter+1;
    model.beta(counter).betaval=array2table(zeros(Init.NbootCV,sum(varaiables_in(j,:))),...
      'VariableNames',Init.FeaturesToOptimizeFrom(varaiables_in(j,:)));
  end
end

disp(['Starting algorithm' model.algorithm])

model.itr=1;
model.state='initialise';
central_message = constructJSONMessage(model,setdiff(fieldnames(model),{'itr'}));

[model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% "Step 2": Optimise models for selection

disp('********Optimiseing models for selection *************')

iter=1;
model.convergence=false;
maxiter = 20;
model.state = 'updateModel';
while ~model.convergence && iter < maxiter
  central_message = constructJSONMessage(model,setdiff(fieldnames(model),{'itr','client'}));
  
  [model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model);
  model = central_update_stratified_cox(model);
  iter=iter+1;
end

%Save result

Result.Init=Init;
Result.Info=model;
temp=[Init.FileSaveParamSelect,'_',datestr(datetime)];
temp=strrep(temp,':','_');
temp=strrep(temp,' ','_');
save(temp,'Result')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Select model using log likelihood

NModels=size(model.client(1).I,1);
NBoots=size(model.client(1).I(1).Result,1);
CVloglik_model = zeros(NBoots,NModels);
for i=1:NModels
  for j=1:length(model.client)
    CVloglik_model(:,i) = CVloglik_model(:,i)+[model.client(j).I(i).Result.CVloglik]';
  end
end
CVLLmean = mean(CVloglik_model);
CVLLstd = std(CVloglik_model);

figure;
%NrParam = sum(model.varselect,2)';
%Calc number of param but correct for levels such that they do not count as
%individual variables
NModelse=size(model.beta,2);
NrParam=zeros(1,NModels);

for i=1:NModelse
  NRemovals=0;
  varnames=model.beta(i).betaval.Properties.VariableNames;
  for j=1:length(Init.FeaturesLevelSets)
    if isempty(setdiff(Init.FeaturesLevelSets{j},varnames))
      NRemovals=NRemovals+1;
      varnames=setdiff(varnames,Init.FeaturesLevelSets{j});
    end
  end
  NrParam(i)=length(varnames)+NRemovals;
end
errorbar(NrParam,CVLLmean,CVLLstd,CVLLstd,'o')
minparam=min(NrParam);
maxparam=max(NrParam);
xlim([minparam-1,maxparam+1]);
indexlow=CVLLmean==min(CVLLmean);
ModelBest=find(indexlow);
ModelBest_varnames=model.beta(ModelBest).betaval.Properties.VariableNames;
cutoff=CVLLmean(indexlow)+CVLLstd(indexlow);
hold on
plot([minparam-1,maxparam+1],[cutoff,cutoff],'r--')
xlabel('Number of parameters in the model')
ylabel('Cross-validated error')
belowcut_index=CVLLmean<=cutoff;
Nrparam1SD=min(NrParam(belowcut_index));
CV1SDmodel=min(CVLLmean(NrParam==Nrparam1SD));
Model1SD=find(CVLLmean==CV1SDmodel & NrParam==Nrparam1SD);
Model1SD_varnames=model.beta(Model1SD).betaval.Properties.VariableNames;

fh1 = figure;
VarName=ModelBest_varnames;
NrParam=length(VarName);
for i=1:NrParam
  subplot(NrParam,1,i)
  histogram(table2array(model.beta(ModelBest).betaval(:,i)))
  title(['Best plot. Param: ',VarName{i}])
end
fh2 = figure;
VarName=Model1SD_varnames;
NrParam=length(VarName);
for i=1:NrParam
  subplot(NrParam,1,i)
  histogram(table2array(model.beta(Model1SD).betaval(:,i)))
  title(['1SD plot. Param: ',VarName{i}])
end

temp=[fullfile(Init.FigurSavePath,'ParameterHist1'),'_',datestr(datetime)];
temp=strrep(temp,':','_');
temp=strrep(temp,' ','_');
savefig(fh1,temp)

temp=[fullfile(Init.FigurSavePath,'ParameterHist2'),'_',datestr(datetime)];
temp=strrep(temp,':','_');
temp=strrep(temp,' ','_');

savefig(fh2,temp)

model.betaval_boot = model.beta;

disp('********Finalised optimiseing models for selection *************')

%Now rum the Model1SD with an increased number of bootstrap
model = rmfield(model,'beta');
model.beta(1).betaval=array2table(zeros(Init.NbootModelInfo,length(Model1SD_varnames)),...
  'VariableNames',Model1SD_varnames);

model.state='terminate_algorithm';
central_message = constructJSONMessage(model,setdiff(fieldnames(model),{'itr','client'}));
[model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model);

iter=1;
%model.itr=1;
model.convergence=false;
maxiter = 20;
model.state = 'updateModel';
while ~model.convergence && iter < maxiter
  central_message = constructJSONMessage(model,setdiff(fieldnames(model),{'itr','client'}));
  
  [model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model);
  model = central_update_stratified_cox(model);
  iter=iter+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% "Step 5": Obtain information from the local servers

disp('********Collect information on model across network *************')

model.alpha=Init.alpha;
model.PI_thresholds = Init.PI_thresholds;
model.CalPlotTimePoints=Init.CalPlotTimePoints;
model.NrCalPlotTimePoints=Init.NrCalPlotTimePoints;
%model.SurvPropTime = [7,15];
%model.SurvPropNrBin = 6;
model.state='localmodelinfo';

%% For lambda = 1
%model.beta.betaval.egelmeer_pi(:)=1;
%% end lambda

central_message = constructJSONMessage(model,setdiff(fieldnames(model),{'itr'}));

[model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model);

disp('********Finished collecting information on model across network *************')

%Save result
Result.Init=Init;
Result.Info=model;
temp=[Init.FileSaveInfo,'_',datestr(datetime)];
temp=strrep(temp,':','_');
temp=strrep(temp,' ','_');
save(temp,'Result')

model.state='terminate_algorithm';
central_message = constructJSONMessage(model,[{'algorithm'},{'state'}]);

[model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model);

model.state='terminate_function';
central_message = constructJSONMessage(model,[{'algorithm'},{'state'}]);
[model, model.itr] = Broadcast_Receive_from_MPI(mpi.svc, central_message, mpi.num_clients, model.itr, model);

%%%% Plotting functions for calibration,  c-Harrell, kaplan meier
modelinfo = cell(length(model.client),1);
for i=1:length(model.client)
  modelinfo{i} = model.client(i).ModelInfo;
end
combinedCalibrationCHarrellKaplanMeierPlots(modelinfo,Init.FigurSavePath)

%Mean of the beta value
disp(mean(table2array(model.beta.betaval)))
toc;
end

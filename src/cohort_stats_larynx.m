function [stats] = cohort_stats_larynx(Data, model)
% This function report back to the central server aggregated information
% about the local data. For numeric data, cumulative distribution and
% values like mean medina std are reported. Within the cumulative
% histogram, the levels are set such that each level at least contains
% NrPtPerBin patients. This value is read from the local config file. A
% typical value would be 5, which ensure that no patient-specific
% information is leaked out of the system. For the categorical variables,
% the number of patients for each level is reported. Data are categorised
% as categorical if the number of unique data values is less than the
% number of patients divided by 10 and the number of levels are not above 5

%Define data types
varnames=Data.Properties.VariableNames;
varnames=setdiff(varnames,'patient');

DateType=false(length(varnames),1);
NumType=false(length(varnames),1);
CatType=false(length(varnames),1);

for i=1:length(varnames)
  DateType(i)=isdatetime(Data.(varnames{i}));
  NumType(i)=isnumeric(Data.(varnames{i}));
  if iscell(Data.(varnames{i}))
    CatType(i)=ischar(Data.(varnames{i}){1});
  else
    if islogical(Data.(varnames{i}))
      CatType(i)=true;
      Data.(varnames{i})=num2str(Data.(varnames{i}));
    end
  end
end
%Check if some of the NumType is a categorical variable. We will assume it
%to be categorical if the uniques number of values are less than the number
%of patients divided by 10. In that case change the numbers to strings
for i=1:length(varnames)
  if NumType(i)
    if length(unique(Data.(varnames{i})))<=size(Data,1)/10 && length(unique(Data.(varnames{i})))<=5
      %Data.(varnames{i})=num2str(Data.(varnames{i}));
      Data.(varnames{i})=arrayfun(@(x) num2str(x),Data.(varnames{i}),'UniformOutput',false);
      NumType(i)=false;
      CatType(i)=true;
    end
  end
end
%Convert all dates to year and treat them as numeric
for i=1:length(varnames)
  if DateType(i)
    Data.(varnames{i})=year(Data.(varnames{i}));
    NumType(i)=true;
    DateType(i)=false;
  end
end
stats.N = size(Data,1);

%Report for the continues data
continuous_var = varnames(NumType);
for i=1:length(continuous_var)
  if any(ismember(Data.Properties.VariableNames,continuous_var{i}))
    stats.(continuous_var{i}).mean = mean(Data.(continuous_var{i})(~isnan(Data.(continuous_var{i}))));
    stats.(continuous_var{i}).median = median(Data.(continuous_var{i})(~isnan(Data.(continuous_var{i}))));
    stats.(continuous_var{i}).std = std(Data.(continuous_var{i})(~isnan(Data.(continuous_var{i}))));
    stats.(continuous_var{i}).sum = sum(Data.(continuous_var{i})(~isnan(Data.(continuous_var{i}))));
    stats.(continuous_var{i}).N = sum(~isnan(Data.(continuous_var{i})));
  end
end
% histograms

NrPtPerBin=model.NrPtPerBin;
for i=1:length(continuous_var)
  NumSamples = sum(~isnan(Data.(continuous_var{i})));
  CumHistLevels=NrPtPerBin/NumSamples;
  CumHistLevels=0:CumHistLevels:1;
  CumHistLevels(end)=1;
  CumHistLevels=prctile(Data.(continuous_var{i})(~isnan(Data.(continuous_var{i}))),100*CumHistLevels);
  CumHist=arrayfun(@(x) sum(Data.(continuous_var{i})(~isnan(Data.(continuous_var{i})))<=x),CumHistLevels);
  %CumHist(1)=0;
  %CumHist(end)=NumSamples;
  CumHist=array2table([CumHistLevels',CumHist'],'VariableNames',{continuous_var{i},'CumulativeNumber'});
  stats.(continuous_var{i}).CumHist=CumHist;
  %stats.(continuous_var{i}).CumHistLevels = CumHistLevels;
end

categorical_variables=varnames(CatType);
for i=1:length(categorical_variables)
  try
    catlevels=unique(Data.(categorical_variables{i}));
  catch
    disp('Stoped in catch statement')
    disp('The stop is at the unique line 87');
    disp(['The variable name is: ',categorical_variables{i}]); 
    if iscell(Data.(categorical_variables{i}))
      disp('The variable class is cell')
      temp=cellfun(@(x) class(x),Data.(categorical_variables{i}),'UniformOutput',false);
      disp('The unique class elements are:');
      disp(unique(temp));
      temp=Data.(categorical_variables{i});
      index=cellfun(@(x) isnumeric(x),Data.(categorical_variables{i}));
      if sum(index)>0
        temp=temp(index);
        temp=unique(temp);
        disp('Unique numeric values:')
        disp(temp)
      end
      index=cellfun(@(x) ischar(x),Data.(categorical_variables{i}));
      if sum(index)>0
        temp=temp(index);
        temp=unique(temp);
        disp('Unique string values:')
        disp(temp)
      end
    else
      disp(['The variable class is: ',class(Data.(categorical_variables{i}))]);
    end
    disp(['The identifier was: ',e.identifier]);
    disp(['The error message was: ',e.message]);
    exit;
  end
  if ~iscell(catlevels)
    catlevels=cellstr(catlevels);
  end
  CategoricalDistibution=cell(size(catlevels,1),2);
  CategoricalDistibution(:,1)=catlevels;
  for j=1:length(catlevels)
    CategoricalDistibution{j,2}=sum(strcmp(Data.(categorical_variables{i}),catlevels(j)));
  end
  CategoricalDistibution=cell2table(CategoricalDistibution,'VariableNames',{categorical_variables{i},'Numbers'});
  stats.(categorical_variables{i})=CategoricalDistibution;
end
%Also, report the number of patients at risk at time points provided from
%the server
PtAtRisk=zeros(length(model.NumerAtRiskTimePoints),1);
for i=1:length(model.NumerAtRiskTimePoints)
  PtAtRisk(i)=sum(Data.survival>=model.NumerAtRiskTimePoints(i)); 
end
PtAtRisk=table(model.NumerAtRiskTimePoints,PtAtRisk,'VariableNames',{'Time','NrPts'});
stats.PtAtRisk=PtAtRisk;
stats.NrPtPerBin=model.NrPtPerBin;
end
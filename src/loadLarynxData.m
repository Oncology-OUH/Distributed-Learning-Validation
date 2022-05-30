function [data, columnHeaders, model]=loadLarynxData(config, model)

data=[]; columnHeaders=[];
datetime.setDefaultFormats('default','yyyy-MM-dd')
if strcmp(config.dataType,'csv')
  
  Data = readtable(config.dataLocation,'Delimiter',',');
  columnHeaders = Data.Properties.VariableNames;
  model.datasourcetype='csv';
  %Data.death_date = datetime(Data.death_date,'InputFormat','dd-MM-yyyy');
  %Data.first_RT_date = datetime(Data.first_RT_date,'InputFormat','dd-MM-yyyy');
  %Data.censor_date = datetime(Data.censor_date,'InputFormat','dd-MM-yyyy');
  
elseif strcmp(config.dataType,'sparql')
  %Config expected to contain endpoint detail for the query sent.
  if isfield(model, 'sparql_query')
    sparql_query = model.sparql_query;
  end
  %[columnHeaders, Data, ~] = sparql_execute(config.endpoint_location, sparql_query, 1);
  [columnHeaders, Data, ~] = sparql_execute(config.dataLocation, sparql_query, 1);
  model.datasourcetype='rdf';
end
%Correct all variable names in data to lower case
Data.Properties.VariableNames=lower(Data.Properties.VariableNames);

%Check all data format and ensure that it is the same format in all
%entries. The main issue is to use a missing format that aligns with the
%format of the variable
DataCell=table2cell(Data);
%Find type for all entrances
temptype=cellfun(@(x) class(x),DataCell,'UniformOutput',false);
DateValues=cellfun(@(x) strcmp(x,'datetime'),temptype);
StringValues=cellfun(@(x) strcmp(x,'char'),temptype);
NumValues=cellfun(@(x) strcmp(x,'double'),temptype);
%Find all missing values
MissingDateValues=cellfun(@(x,y) y && isnan(x.Year),DataCell,num2cell(DateValues));
MissingStringValues=cellfun(@(x,y) y && (strcmp(x,'') || strcmp(x,' ') || strcmp(x,'.') || strcmp(x,'-')),DataCell,num2cell(StringValues));
MissingNumValues=cellfun(@(x,y) isnumeric(x) && (logical(sum(isnan(x))) || isempty(x)),DataCell,num2cell(NumValues));
DataCell(MissingStringValues)={''};

MissingValues=MissingDateValues | MissingStringValues | MissingNumValues;

datecollumns=repmat(sum(DateValues)>0,size(DataCell,1),1);
numcollumns=repmat(sum(StringValues & ~MissingValues)==0,size(DataCell,1),1) & ~datecollumns;
stringcollumns=repmat(sum(NumValues & ~MissingValues)==0,size(DataCell,1),1) & ~datecollumns;
mixcollumns=~datecollumns & ~numcollumns & ~stringcollumns;

%Date columns make all non-dates to NAT
DataCell(datecollumns & ~DateValues)={NaT};
%Num columns make all missing to an empty number
DataCell(numcollumns & MissingValues)={NaN};
%String columns set non-string to the empty string
DataCell(stringcollumns & MissingValues)={''};
%Mix columns, change numbers to string and set missing to empty string
DataCell(mixcollumns & NumValues)=cellfun(@(x) num2str(x),DataCell(mixcollumns & NumValues),'UniformOutput',false);
DataCell(mixcollumns & MissingValues)={''};

Data=cell2table(DataCell,'VariableNames',Data.Properties.VariableNames);
%End check all data format

%Check and correct format of dates
DateVariables={'death_date','first_rt_date','censor_date'};
for j=1:length(DateVariables)
  for i=1:length(Data.(DateVariables{j}))
    if ~isdatetime(Data.(DateVariables{j})(i))
      if iscell(Data.(DateVariables{j})(i))
        if ~isdatetime(Data.(DateVariables{j}){i})
          if ischar(Data.(DateVariables{j}){i})
            temp=strfind(Data.(DateVariables{j}){i},'-');
            if length(temp)==2 && temp(1)==3 && temp(2)==6
              if strcmp(Data.(DateVariables{j})(i),'NA')
                Data.(DateVariables{j}){i}=NaT;
              else
                Data.(DateVariables{j}){i}=datetime(Data.(DateVariables{j})(i),'InputFormat','dd-MM-yyyy');
              end
            else
              Data.(DateVariables{j}){i}=datetime(Data.(DateVariables{j})(i),'InputFormat','yyyy-MM-dd');
            end
          else
            Data.(DateVariables{j}){i}=NaT;
          end
        end
      else
        Data.(DateVariables{j}){i}=NaT;
      end
    end
  end
  if iscell(Data.(DateVariables{j}))
    Data.(DateVariables{j})=table2array(cell2table(Data.(DateVariables{j})));
  end
end
Data.vital_status = ~isnat(Data.death_date);

%survival time in months
for i=1:length(Data.vital_status)
  if ~Data.vital_status(i)
    Data.survival(i) = ((hours(Data.censor_date(i) - Data.first_rt_date(i))/24)/365.25)*12;
  else
    Data.survival(i) = ((hours(Data.death_date(i) - Data.first_rt_date(i))/24)/365.25)*12;
  end
end
if ~isempty(model.HardCensoringTime)
  index=Data.survival>model.HardCensoringTime;
  Data.survival(index)=model.HardCensoringTime;
  Data.vital_status(index)=0;
end
%Remove all patients with T stage Tis which is not cancer - yet
if isnumeric(Data.t_stage)
  index=contains(Data.t_stage,'Tis','IgnoreCase',true);
  Data=Data(~index,:);
end
%Remove all patients with first RT date outside the years 2005 to 2018
index=Data.first_rt_date.Year>=2005 & Data.first_rt_date.Year<=2018;
Data=Data(index,:);

%Remove data outside the requested dose and fraction range (only if the data is available in the data)
if sum(strcmp('dose',Data.Properties.VariableNames))>=1
  if isfield(model,'TotalDoseRange')
    index=Data.dose>=model.TotalDoseRange(1) & Data.dose<=model.TotalDoseRange(2);
  end
end
Data=Data(index,:);
if sum(strcmp('fractions',Data.Properties.VariableNames))>=1
  if isfield(model,'TotalFracRange')
    index=Data.fractions>=model.TotalFracRange(1) & Data.fractions<=model.TotalFracRange(2);
  end
end
Data=Data(index,:);


model.stats_all = cohort_stats_larynx(Data, model);

%Reduce the dataset to include only those patient that has at less than or
%equal to the number of allowed missing values.
temp=table2cell(Data(:,model.ImputationVariables));
MissingValues=cellfun(@(x) logical(sum(isnan(x))) || isempty(x),temp);
N_missing=sum(MissingValues,2);
index=N_missing<=model.NAllowedMissing;
model.stats_AllowedMissing = cohort_stats_larynx(Data(index,:), model);
model.stats_NoMissing = cohort_stats_larynx(Data(N_missing==0,:), model);
%Change T sublevels to the four main levels
%index1=contains(Data.t_stage,'1');
if isnumeric(Data.t_stage)
  Data.t_stage=num2str(Data.t_stage);
end
index2=contains(Data.t_stage,'2');
index3=contains(Data.t_stage,'3');
index4=contains(Data.t_stage,'4');
index5=contains(Data.t_stage,'TX','IgnoreCase',true);
index1=~(index2 | index3 | index4 | index5 |strcmp(Data.t_stage,''));
Data.t_stage(index1)={'T1'};
Data.t_stage(index2)={'T2'};
Data.t_stage(index3)={'T3'};
Data.t_stage(index4)={'T4'};
Data.t_stage(index5)={''};

%Merge N1 N2, and N3 since we are only interested in N0 versus N+
if isnumeric(Data.n_stage)
  Data.n_stage=num2str(Data.n_stage);
end
%index0=contains(Data.n_stage,'0');
index1=contains(Data.n_stage,'1');
index2=contains(Data.n_stage,'2');
index3=contains(Data.n_stage,'3');
index4=contains(Data.n_stage,'NX','IgnoreCase',true);
index0=~(index1 | index2 | index3 | index4 | strcmp(Data.n_stage,''));
Data.n_stage(index0)={'N0'};
Data.n_stage(~index0 & ~strcmp(Data.n_stage,''))={'N1'};
Data.n_stage(index4)={''};

indexGlottis = contains(Data.tumour_loc,'Glottis Carcinoma');
Data.tumour_loc(indexGlottis)={'Glottis Carcinoma'};
Data.tumour_loc(~indexGlottis & ~strcmp(Data.tumour_loc,''))={'Non-Glottis'};

Data=Data(index,:);
MissingValues=sortrows(MissingValues);
model.MissingValues=array2table(MissingValues,'VariableNames',model.ImputationVariables);

NBoot=size(model.beta(1).betaval,1);
NPts=size(Data,1);
rng(model.randseed+model.localseed);
randomvalues=ceil(NPts*rand(NPts,NBoot));
DataImputed=SimpleMice_OUH('Data',Data(:,model.ImputationVariables),'DataTypes',struct2table(model.ImputationVariablesTypes),'NrDataSets',NBoot,'Seed',model.randseed+model.localseed);

ImputedDataTrain=cell(size(randomvalues,2),1);
ImputedDataValid=cell(size(randomvalues,2),1);
for i=1:NBoot
  temp=Data;
  temp(:,model.ImputationVariables)=DataImputed{i}(:,model.ImputationVariables);
  
  temp.genderMale=ones(size(Data,1),1);
  temp.genderMale(strcmpi(temp.gender,'Female'))=0;
  
  temp.T2 = double(contains(temp.t_stage,'2'));
  temp.T3 = double(contains(temp.t_stage,'3'));
  temp.T4 = double(contains(temp.t_stage,'4'));
  temp.Nplus = double(~contains(temp.n_stage,'0'));
  temp.NonGlottic = double(contains(temp.tumour_loc,'Non-Glottis'));
  
  egelmeer_data = temp(:,fieldnames(model.PI_coefficients));
  %Data.egelmeer_pi = table2array(egelmeer_data)*model.egelmeer_beta2yr(:);
  temp.egelmeer_pi = table2array(egelmeer_data)*(struct2array(model.PI_coefficients)');
  
  
  indextrain=randomvalues(:,i);
  indexvalid=setdiff(1:NPts,indextrain);
  ImputedDataTrain{i}=temp(indextrain,:);
  ImputedDataValid{i}=temp(indexvalid,:);
  
end
data.ImputedDataTrain=ImputedDataTrain;
data.ImputedDataValid=ImputedDataValid;
end

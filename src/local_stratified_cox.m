function [local_message,model,mpi,itr]=local_stratified_cox(mpi,itr,central_message,model,data)

local_message=[];
%I do not see that the next two lines are used, and model.state will be
%overwritten by the json call just after these lines. I have added and
%initialisation just before the json call to ensure that it is always
%declared independent of the result of the json call
%message_tokens = strsplit(central_message,';');
%model.state = message_tokens(ismember(message_tokens,'state')+1);
model.state='';
[model, ~] = readJSONMessage('local',model,central_message,model.alg_client_index);



if strcmp(model.state,'updateModel')
  
  
  %if model.NModel==1
  %     model.varselect = model.varselect';
  %end
  
  Nboot = size(model.beta(1).betaval,1);
  NModel = size(model.beta,1);
  %   for j=1:NModel
  %     model.loglik2(j).L = zeros(Nboot,1);
  %     model.CVloglik2(j).L = zeros(Nboot,1);
  %   end
  %model.loglik = zeros(NModel*Nboot,1);
  %model.CVloglik = zeros(NModel*Nboot,1);
  %model.gradient = zeros(NModel*Nboot,size(model.varselect,2));
  %model.hessian = zeros(NModel*Nboot,size(model.varselect,2)^2);
  centerval=cell(NModel,1);
  for j=1:NModel
    
    betaval=struct2table(model.beta(j).betaval);
    betavarnames=betaval.Properties.VariableNames;
    I=cell(Nboot,1);
    %The following line can be changed to a parfor line if executed on a multi-core machine locally to increase the calculation speed of the bootstrap method
    for i=1:Nboot
      xvalTemp=data.ImputedDataTrain{i}(:,(betavarnames));
      
      SurvivalTemp=data.ImputedDataTrain{i}(:,{'survival','vital_status'});
      SurvivalTemp.Properties.VariableNames={'Time','Event'};
      
      I{i}=PartialLikelihoodCox(SurvivalTemp,xvalTemp,betaval(i,:));
      temp=PartialLikelihoodCox(SurvivalTemp,xvalTemp,array2table(zeros(size(betaval(i,:))),'VariableNames',betaval.Properties.VariableNames));
      
      I{i}.loglik=temp.loglik-I{i}.loglik;
      
      %model.loglik(((j-1)*Nboot)+i) = temp.loglik - I{i}.loglik;
      %model.gradient(((j-1)*Nboot)+i,:) = [I{i}.gradient(:)' zeros(1,size(model.gradient,2)-length(I{i}.gradient(:)'))];
      %model.hessian(((j-1)*Nboot)+i,:) = [I{i}.hessian(:)' zeros(1,size(model.hessian,2)-length(I{i}.hessian(:)'))];
      %Check that there is data in the validation cohort
      if sum(size(data.ImputedDataValid{i},1))>0
        xvalTemp=data.ImputedDataValid{i}(:,(betavarnames));
        SurvivalTemp=data.ImputedDataValid{i}(:,{'survival','vital_status'});
        SurvivalTemp.Properties.VariableNames={'Time','Event'};
        
        temp1=PartialLikelihoodCox(SurvivalTemp,xvalTemp,betaval(i,:));
        temp2=PartialLikelihoodCox(SurvivalTemp,xvalTemp,array2table(zeros(size(betaval(i,:))),'VariableNames',betaval.Properties.VariableNames));
        
        %A reference to the scaling by the number of patients in the next
        %line (not the number of events) is according to F. Harrell “Regression Modeling Strategies” (page 505 REF 548). This reference
        %is a three-page article by M. SCHEMPER Biometrika, Volume 79, Issue 1, March 1992, Pages 202–204.
        I{i}.CVloglik=(temp2.loglik-temp1.loglik)*size(data.ImputedDataTrain{i},1)/size(xvalTemp,1);
        %model.CVloglik(((j-1)*Nboot)+i) = (temp2.loglik-temp1.loglik)*size(xval,1)/size(xvalTemp,1);
      else
        I{i}.CVloglik=NaN;
        %model.CVloglik(((j-1)*Nboot)+i) = NaN;
      end
      
    end
    temp=[];
    temp.VarNames=betavarnames;
    temp.Result=I;
    centerval{j}=temp;
  end
  model.I=centerval;
  local_message = {'client_name','I'};
  
  
  
elseif strcmp(model.state,'localmodelinfo')
  try
    %Within the localmodelinfo we will only report data on the first model if
    %more models are forwarded for analysis (of course, all bootstraps for
    %that model are included). This is to prevent “the harvesting of all
    %data” in just one run.
    model.beta.betaval=struct2table(model.beta(1).betaval);
    
    betavarnames=model.beta.betaval.Properties.VariableNames;
    
    NBoot=size(model.beta.betaval,1);
    
    betaval=table2array(model.beta.betaval)';
    %Calculate the distribution of PI individually for each set of imputed
    %datasets and related beta values.
    LinPredict_train=cell(NBoot,1);
    LinPredict_train_medianbeta=cell(NBoot,1);
    LinPredict_valid=cell(NBoot,1);
    LinPredict_valid_medianbeta=cell(NBoot,1);
    
    %Calculate the median beta. The median beta can e.g. be used to calculate
    %the uncertainty of the baseline hazard. The uncertainty based on the
    %median value is the uncertainty that is of interest in case baseline
    %differences between institutions are of interest. In contrast, the
    %individual beta will provide the overall uncertainty of the baseline
    %hazard, which then will include an uncertainty related to the
    %uncertainty in beta. The reason for this is that a high beta will lead
    %to a low baseline and vice versa; thus, they are directly linked to each
    %other
    
    MedianBeta=median(betaval,2);
    for i=1:NBoot
      xval=data.ImputedDataTrain{i}(:,betavarnames);
      temp=table2array(xval);
      LinPredict_train{i}=temp*(betaval(:,i));
      LinPredict_train_medianbeta{i}=temp*MedianBeta;
      xval=data.ImputedDataValid{i}(:,betavarnames);
      temp=table2array(xval);
      LinPredict_valid{i}=temp*(betaval(:,i));
      LinPredict_valid_medianbeta{i}=temp*MedianBeta;
    end
    
    
    %Start histogram of PI
    LinPreVersions_train={LinPredict_train,LinPredict_train_medianbeta};
    LinPreLabel={'IndividualBeta','MedianBeta'};
    for j=1:length(LinPreVersions_train)
      temp=LinPreVersions_train{j};
      AllLinPreTrain = cat(1,temp{:});
      CumHistLevels=NBoot*model.NrPtPerBin/size(AllLinPreTrain,1);
      CumHistLevels=min(CumHistLevels,1);
      CumHistLevels=0:CumHistLevels:1;
      CumHistLevels(end)=1;
      CumHistLevels=prctile(AllLinPreTrain,100*CumHistLevels);
      CumHist=zeros(NBoot,length(CumHistLevels));
      for i=1:NBoot
        CumHist(i,:)=arrayfun(@(x) sum(temp{i}<=x),CumHistLevels);
      end
      
      temp=prctile(CumHist,100*[(model.alpha/2),(1-model.alpha/2),.5],1);
      PI=table(CumHistLevels',temp(3,:)',temp(1,:)',temp(2,:)','VariableNames',...
        {'CumHist_Levels','CumHist','CumHist_Lower','CumHist_Upper'});
      I.PI.(LinPreLabel{j})=PI;
    end
    %End histogram of model variables and PI
    
    %Start calculating the baseline hazard and C-Harrell index
    
    
    %Start calculating confidence interval baseline hazard and C_Harrell index
    LinPreVersions_train={LinPredict_train,LinPredict_train_medianbeta};
    LinPreVersions_valid={LinPredict_valid,LinPredict_valid_medianbeta};
    LinPreLabel={'IndividualBeta','MedianBeta'};
    for z=1:length(LinPreVersions_train)
      temp=cellfun(@(x) x.survival,data.ImputedDataTrain,'UniformOutput',false);
      UniqueTimes = unique(cat(1,temp{:}));
      HazardBootRes=zeros(length(UniqueTimes),NBoot);
      CHarrellBootRes=zeros(NBoot,1);
      CHarrellBootOutRes=zeros(NBoot,1);
      for i=1:NBoot
        %disp(['z: ',num2str(z),' i: ',num2str(i)])
        %LinPredict_boot=LinPredict_train{i};
        LinPredict_boot=LinPreVersions_train{z}{i};
        tempptdata_inboot=table2array(data.ImputedDataTrain{i}(:,{'survival','vital_status'}));
        tempptdata_inboot=array2table([tempptdata_inboot,LinPredict_boot],'VariableNames',{'Time','Event','LinPredict'});
        
        Hazard_boot_temp=CalcHazard(tempptdata_inboot);
        HazardBootRes(:,i)=interp1(Hazard_boot_temp.Time,Hazard_boot_temp.CumHazard,UniqueTimes,'next');
        
        CHarrellBootRes(i)=CHarrellIndex(tempptdata_inboot.Time,tempptdata_inboot.Event,LinPredict_boot);
        LinPredict_outboot=LinPreVersions_valid{z}{i};
        tempptdata_outboot=table2array(data.ImputedDataValid{i}(:,{'survival','vital_status'}));
        tempptdata_outboot=array2table([tempptdata_outboot,LinPredict_outboot],'VariableNames',{'Time','Event','LinPredict'});
        CHarrellBootOutRes(i)=CHarrellIndex(tempptdata_outboot.Time,tempptdata_outboot.Event,LinPredict_outboot);
        
      end
      
      temp=prctile(CHarrellBootRes,100*[(model.alpha/2),(1-model.alpha/2),.5]);
      CHarrell.CHarrellBootLower=temp(1);
      CHarrell.CHarrellBootUpper=temp(2);
      CHarrell.CHarrellBootMedian=temp(3);
      CHarrell.Optimism=mean(CHarrellBootRes-CHarrellBootOutRes);
      temp=prctile(CHarrellBootOutRes,100*[(model.alpha/2),(1-model.alpha/2),.5]);
      CHarrell.CV_CHarrellBootLower=temp(1);
      CHarrell.CV_CHarrellBootUpper=temp(2);
      CHarrell.CV_CHarrellBootMedian=temp(3);
      
      I.CHarrell.(LinPreLabel{z})=CHarrell;
      
      
      
      Confidence=zeros(size(HazardBootRes,1),3);
      for i=1:size(HazardBootRes,1)
        Confidence(i,:)=prctile(HazardBootRes(i,:),100*[(model.alpha/2),(1-model.alpha/2),.5]);
      end
      tempname = Hazard_boot_temp.Properties.VariableNames;
      
      Hazard_model=table(UniqueTimes,Confidence(:,3),Confidence(:,1),Confidence(:,2));
      Hazard_model.Properties.VariableNames = [tempname,{'Lower'},{'Upper'}];
      
      temp=cellfun(@(x) x.survival,data.ImputedDataTrain,'UniformOutput',false);
      AllSurivivalTimes=cat(1,temp{:});
      TimeBins=NBoot*model.NrPtPerBin/size(AllSurivivalTimes,1);
      TimeBins=min(TimeBins,1);
      TimeBins=0:TimeBins:1;
      TimeBins(end)=1;
      TimeBins=prctile(AllSurivivalTimes,100*TimeBins);
      TimeBins=TimeBins';
      
      I.BaselineHazard.(LinPreLabel{z})=[];
      I.BaselineHazard.(LinPreLabel{z}).Time=TimeBins;
      I.BaselineHazard.(LinPreLabel{z}).CumHazard=interp1(Hazard_model.Time,Hazard_model.CumHazard,TimeBins,'next');
      I.BaselineHazard.(LinPreLabel{z}).Lower=interp1(Hazard_model.Time,Hazard_model.Lower,TimeBins,'next');
      I.BaselineHazard.(LinPreLabel{z}).Upper=interp1(Hazard_model.Time,Hazard_model.Upper,TimeBins,'next');
      I.BaselineHazard.(LinPreLabel{z})=struct2table(I.BaselineHazard.(LinPreLabel{z}));
      
      %End calculating the baseline hazard and C-Harrell index
      
      %%%% Kaplan Meier
      if isstruct(model.PI_thresholds)
        %model.PI_thresholds has center individual cutvalues
        centername=fieldnames(model.PI_thresholds);
        index=strcmp(model.client_name,centername);
        if sum(index)==0
          index(1)=true;
        end
        centername=centername(index);
        local_PI_thresholds=model.PI_thresholds.(centername{1});
      else
        %Same cutvalues for all centers
        local_PI_thresholds=model.PI_thresholds;
      end
      local_PI_thresholds=unique([-Inf, local_PI_thresholds(:)',Inf]);
      %model.PI_thresholds=unique([-Inf, model.PI_thresholds(:)',Inf]);
      
      KM_Cox.(LinPreLabel{z})=cell(size(local_PI_thresholds,2)-1,1);
      for i=1:(size(local_PI_thresholds,2)-1)
        KM=zeros(length(UniqueTimes),NBoot);
        Events=zeros(length(UniqueTimes),NBoot);
        Censors=zeros(length(UniqueTimes),NBoot);
        CoxSurvProp=zeros(length(UniqueTimes),NBoot);
        for j=1:NBoot
          %disp(['i: ',num2str(i),' j: ',num2str(j)])
          LinPre=LinPreVersions_train{z}{j};
          %LinPre=LinPredict_train{j};
          index_PI_bin=(LinPre>=local_PI_thresholds(i)) & (LinPre<local_PI_thresholds(i+1));
          ptdataPIBin=data.ImputedDataTrain{j}(index_PI_bin,{'survival','vital_status'});
          index=ptdataPIBin.vital_status;
          if sum(index_PI_bin)>=2 && size(unique(ptdataPIBin(index,'survival')),1)>1
            %There is at least two events with different event times
            ptdataPIBin.Properties.VariableNames={'Time','Event'};
            [ftemp,x,~,~]=KaplanMeierEstimator(ptdataPIBin,model.alpha);
            
            KM(:,j)=interp1(x,ftemp,UniqueTimes,'previous');
            %KM(x<min(UniqueTimes))=1;
            %KM(x>max(UniqueTimes))=0;
            KM(UniqueTimes<min(x),j)=1;
            if max(UniqueTimes)>max(x)
              %For bootstraps for which the last time point is less than the
              %maximum time point within all the data, make an extrapolation
              %of the survival curve, using a Weibull distribution, to get an
              %unbiased evaluation across all bootstrapped values
              indexpos=ptdataPIBin.Time>0;
              Weibullfit = fitdist(ptdataPIBin.Time(indexpos),'wbl','Censoring',~ptdataPIBin.Event(indexpos));
              S=1-cdf('wbl',UniqueTimes,Weibullfit.A,Weibullfit.B);
              KM(UniqueTimes>max(x),j)=min(ftemp(end),S(UniqueTimes>max(x)));
            end
            %KM(UniqueTimes>max(x),j)=ftemp(end);
            for k=1:length(UniqueTimes)
              if k==1
                index=(ptdataPIBin.Time<=UniqueTimes(k));
              else
                index=(ptdataPIBin.Time>UniqueTimes(k-1)) & (ptdataPIBin.Time<=UniqueTimes(k));
              end
              Events(k,j)=sum(ptdataPIBin.Event(index));
              Censors(k,j)=sum(index)-Events(k,j);
            end
            Events(:,j)=cumsum(Events(:,j),1);
            Censors(:,j)=cumsum(Censors(:,j),1);
            for k=1:length(UniqueTimes)
              CumHazard_0=interp1(UniqueTimes,HazardBootRes(:,j),UniqueTimes(k),'next');
              CoxSurvProp(k,j)=mean(exp(-CumHazard_0*exp(LinPre(index_PI_bin))));
            end
            
          else
            KM(:,j)=NaN;
            Events(:,j)=NaN;
            Censors(:,j)=NaN;
            CoxSurvProp(:,j)=NaN;
          end
        end %Boot
        
        temp=prctile(KM,100*[(model.alpha/2),(1-model.alpha/2),.5],2);
        KM_Cox.(LinPreLabel{z}){i}.LowerPI=local_PI_thresholds(i);
        KM_Cox.(LinPreLabel{z}){i}.UpperPI=local_PI_thresholds(i+1);
        KM_Cox.(LinPreLabel{z}){i}.Values.Time = TimeBins;
        KM_Cox.(LinPreLabel{z}){i}.Values.KM = interp1(UniqueTimes,temp(:,3),TimeBins,'previous');
        KM_Cox.(LinPreLabel{z}){i}.Values.KM_Lower = interp1(UniqueTimes,temp(:,1),TimeBins,'previous');
        KM_Cox.(LinPreLabel{z}){i}.Values.KM_Upper = interp1(UniqueTimes,temp(:,2),TimeBins,'previous');
        temp=prctile(Events,100*[(model.alpha/2),(1-model.alpha/2),.5],2);
        KM_Cox.(LinPreLabel{z}){i}.Values.CumEvents = interp1(UniqueTimes,temp(:,3),TimeBins,'previous');
        KM_Cox.(LinPreLabel{z}){i}.Values.CumEvents_Lower = interp1(UniqueTimes,temp(:,1),TimeBins,'previous');
        KM_Cox.(LinPreLabel{z}){i}.Values.CumEvents_Upper = interp1(UniqueTimes,temp(:,2),TimeBins,'previous');
        temp=prctile(Censors,100*[(model.alpha/2),(1-model.alpha/2),.5],2);
        KM_Cox.(LinPreLabel{z}){i}.Values.CumCensoring = interp1(UniqueTimes,temp(:,3),TimeBins,'previous');
        KM_Cox.(LinPreLabel{z}){i}.Values.CumCensoring_Lower = interp1(UniqueTimes,temp(:,1),TimeBins,'previous');
        KM_Cox.(LinPreLabel{z}){i}.Values.CumCensoring_Upper = interp1(UniqueTimes,temp(:,2),TimeBins,'previous');
        temp=prctile(CoxSurvProp,100*[(model.alpha/2),(1-model.alpha/2),.5],2);
        KM_Cox.(LinPreLabel{z}){i}.Values.Cox = interp1(UniqueTimes,temp(:,3),TimeBins,'previous');
        KM_Cox.(LinPreLabel{z}){i}.Values.Cox_Lower = interp1(UniqueTimes,temp(:,1),TimeBins,'previous');
        KM_Cox.(LinPreLabel{z}){i}.Values.Cox_Upper = interp1(UniqueTimes,temp(:,2),TimeBins,'previous');
        KM_Cox.(LinPreLabel{z}){i}.Values=struct2table(KM_Cox.(LinPreLabel{z}){i}.Values);
        
      end %PI
      
      
      I.KM_Cox=KM_Cox;
      %Calibration plots
      NPatients=size(data.ImputedDataTrain{1},1);
      
      for i=1:length(model.CalPlotTimePoints)
        x=[];
        x.time=model.CalPlotTimePoints(i);
        PredictedSurvival=zeros(NPatients,NBoot);
        for j=1:NBoot
          CumHazard_0=interp1(UniqueTimes,HazardBootRes(:,j),model.CalPlotTimePoints(i),'next');
          PredictedSurvival(:,j)=exp(-CumHazard_0*exp(LinPreVersions_train{z}{j}));
        end
        temp=0:1/model.NrCalPlotTimePoints:1; temp(end)=1;
        SurvivalCutPoints=prctile(PredictedSurvival(:),100*temp);
        
        BinNumber = interp1(SurvivalCutPoints,1:(model.NrCalPlotTimePoints+1),PredictedSurvival,'previous');
        
        BinNumber(BinNumber==model.NrCalPlotTimePoints+1)=model.NrCalPlotTimePoints;
        KMtemp=zeros(model.NrCalPlotTimePoints,NBoot);
        for j=1:NBoot
          ptdata=data.ImputedDataTrain{j}(:,{'survival','vital_status'});
          ptdata.Properties.VariableNames={'Time','Event'};
          for k=1:model.NrCalPlotTimePoints
            index=BinNumber(:,j)==k;
            if sum(index)>=1 && sum(ptdata(index,:).Event)>0
              [f,timeval,~,~]=KaplanMeierEstimator(ptdata(index,:),model.alpha);
              KMtemp(k,j)=interp1(timeval,f,model.CalPlotTimePoints(i),'previous');
              if model.CalPlotTimePoints(i)<min(timeval)
                KMtemp(k,j)=1;
              end
              if model.CalPlotTimePoints(i)>max(timeval)
                %For bootstraps for which the last time point is less than the
                %requested time point, make an extrapolation
                %of the survival curve, using a Weibull distribution, to get an
                %unbiased evaluation across all bootstrapped values
                indexpos=ptdata.Time>0;
                Weibullfit = fitdist(ptdata(index&indexpos,:).Time,'wbl','Censoring',~ptdata(index&indexpos,:).Event);
                S=1-cdf('wbl',model.CalPlotTimePoints(i),Weibullfit.A,Weibullfit.B);
                KM(k,j)=min(f(end),S);
              end
            else
              KMtemp(k,j)=NaN;
            end
          end
        end
        Predicted=zeros(model.NrCalPlotTimePoints,1);
        
        for k=1:model.NrCalPlotTimePoints
          index=BinNumber==k;
          
          Predicted(k)=median(PredictedSurvival(index(:)));
        end
        x.Plot.Predicted = Predicted;
        temp=prctile(KMtemp,100*[model.alpha/2,1-model.alpha/2,.5],2);
        x.Plot.KM_Observed = temp(:,3);
        x.Plot.Lower = temp(:,1);
        x.Plot.Upper = temp(:,2);
        x.Plot=struct2table(x.Plot);
        I.CalPlot.(LinPreLabel{z}){i}=x;
      end
    end  %EndLinPre version
    
    model.ModelInfo = I;
    local_message = {'client_name','ModelInfo'};
  catch e
    save('model','model');
    save('data','data');
    disp('Stoped in catch statement')
    disp(['The identifier was: ',e.identifier]);
    disp(['The error message was: ',e.message]);
    exit;
  end
elseif strcmp(model.state,'initialise')
  local_message = {'client_name','stats_all','stats_AllowedMissing','stats_NoMissing','MissingValues','datasourcetype'};
else
  disp('Model state not found.')
end


end
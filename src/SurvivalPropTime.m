function I=SurvivalPropTime(ptdata,Hazard_model,timepoint,nrbin,alpha)
I=cell(length(timepoint),1);
for i=1:length(timepoint)
  CumHazard_0=interp1(Hazard_model.Time,Hazard_model.CumHazard,timepoint(i),'previous'); %Baseline hazard at timepoint
  SurvProp=exp(-CumHazard_0*exp(ptdata.LinPredict)); %Predicted survival at timepoint
  SurvivedTime=NaN(size(ptdata,1),1); %Array used to indicate whether the patient survived beyond timepoint
  index=ptdata.Time>timepoint(i);
  SurvivedTime(index)=1;
  index=ptdata.Time<=timepoint(i) & ptdata.Event==1;
  SurvivedTime(index)=0; %After this step, the observed survival after timepoint contain 0, 1, or NaN (NaN if censored before timepoint)
  x=[];
  x.time=timepoint(i);
  indexSurvivalKnown=~isnan(SurvivedTime);
  DistSurvProp=SurvProp(indexSurvivalKnown); %Predicted survival for those patients for which their survival beyond timepoint is known
  temp=0:1/nrbin:1; temp(end)=1;
  CutPointSurvProp=prctile(DistSurvProp,100*temp);
  
  BinNumber = interp1(CutPointSurvProp,1:(nrbin+1),SurvProp,'previous'); %Bin risk number for patients
  calplotdata=zeros(nrbin,4);
  for j=1:nrbin
    index_bin=BinNumber==j;
    [f,timeval,flo,fup]=KaplanMeierEstimator(ptdata(index_bin,:),alpha);

    calplotdata(j,1)=median(SurvProp(index_bin));
    calplotdata(j,2)=interp1(timeval,f,timepoint(i),'previous');
    calplotdata(j,3)=interp1(timeval,flo,timepoint(i),'previous');
    calplotdata(j,4)=interp1(timeval,fup,timepoint(i),'previous');
  end
  x.Plot.Predicted = calplotdata(:,1);
  x.Plot.KM_Observed = calplotdata(:,2);
  x.Plot.Lower = calplotdata(:,3);
  x.Plot.Upper = calplotdata(:,4);
  I{i}=x;
end
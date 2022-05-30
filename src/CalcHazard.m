function I=CalcHazard(data)
%The estimator formula can be found in the book “Regression Modeling Strategies” by Frank E. Harrel Jr on page 485
% The returned value is the cumulative hazard function. Thus, in the
% comments below, when the word hazard is used, it is implicitly understood
% to mean the cumulative hazard.
uniquetime=unique(data.Time);
Hazard=zeros(length(uniquetime),1);
for i=1:length(uniquetime)
  largetime_index=data.Time>=uniquetime(i);
  sumlarge=sum(exp(data.LinPredict(largetime_index)));
  death=sum(data.Event(data.Time==uniquetime(i)));
  if sumlarge>0
    Hazard(i)=death/sumlarge;
  else
    Hazard(i)=0;
  end
end
Hazard=cumsum(Hazard);
%The hazard function is a step function that is constant between the time
%points. At a given event or censoring time, the hazard function has the
%value as just before the event or censoring. The hazard is changed just
%after the timepoint if an event occurred (not if it was censoring). Thus,
%the individual hazard intervals are open toward negative times and closed
%toward positive time. Therefore the hazard at the first event time will be
%zero. If one wants to know the value of the hazard function just after
%that event, it is needed to find the hazard function value at the next
%point in time.

%The open and closed ends of the intervals of the hazard intervals are the
%opposite for the Kaplan-Meier estimator, in which the Survival at the time
%of the first event is less than unity. The idea is that the hazard is the
%same until the event is observed and then changed afterwards. Thus for the
%Kaplan Meier estimator, one needs to look for the previous timepoint to
%find the value just before a specific timepoint. This needs to be kept in
%mind when interpolation KM and Hazards.

Hazard=[0; Hazard(1:end-1)]; %This correction is due to the summation over time in the estimator should be less than the time of interest

%Since the hazard is zero until the first event, it is always okay to add a
%hazard of zero at time zero. It actually only states that hazard up to and
%including that timepoint (all negative times and zero) has a hazard of
%zero. That is not very informative since that is the case by definition.
%However, if one wants to use Matlab to interpolate hazard values, it might
%be helpful to include the timepoint zero explicitly.
if uniquetime(1)>0           
  uniquetime=[0;uniquetime];
  Hazard=[0;Hazard];
end  
I=array2table([uniquetime,Hazard],'VariableName',{'Time','CumHazard'});
end
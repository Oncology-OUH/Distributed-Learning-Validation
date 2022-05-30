function I=CHarrellIndex(Time,Event,LinPredictor)
% In “Evaluating the Yield of Medical Test” by Harrell et al Jama May 14
% 1982 Vol 247 No 18 Harrells C-index for survival model (Cox model using a
% linear predictor (ln(-ln(P))= b1x1+b2x2+…+bnxn = PI = linear predictor))
% is described in the following way: Draw a pair of patients and determine
% which patient lived longer from his baseline evaluation. Survival times
% can be validly compared either when both patients have died or when one
% has died and the other’s follow-up time has exceeded the survival time of
% the first. If both patients are still alive, which will live the longer
% is not known, and that pair of patients is not used in the analysis.
% Otherwise, it can be determined whether the patient with the higher
% prognostic score (e.g. low linear predictor) also had the longer survival
% time. The process is repeated until all possible pairs of patients have
% been examined. Of the pairs of patients for which the ordering of
% survival times could be inferred, the fraction of pairs such that the
% patient with the higher score (e.g. low linear predictor) had the longer
% survival time will be denoted by c.
%
% The current program implements this calculation and return Harrells C-index.
%
% The input parameters are the
% survival times (e.g. follow-up times) in a 1D array, a matching array
% indicating whether the follow-up time is and event (value 1 or true) (0
% or false being censored values), and in the third array the values of the
% linear predictor related to the patients.

% CaB 19/10-2018

I=0;
if size(Time,1)==1
  Time=Time';
end
if size(Event,1)==1
  Event=Event';
end
if size(LinPredictor,1)==1
  LinPredictor=LinPredictor';
end
if ((size(Time,1) ~= size(Event,1)) || (size(Event,1)~=size(LinPredictor,1)))
  %The three input arrays differ in length
  warning('The input arrays do not have the same length which is needed for this calculation');
  return;
end
if size(Time,2)~=1 || size(Event,2)~=1 || size(LinPredictor,2)~=1
  warning('One of the three input arrays are not 1D arrays which is needed for this calculation');
  return;
end

%Convert Event to logical if it is not delivered as logical
temp=Event;
if ~islogical(temp)
  Event=false(size(Event,1),1);
  Event(temp==1)=true;
end
if sum(Event)==0
  warning('There is no events which is needed to perform the calculation of the Harrell index');
  return;
end


Numerator=0;
Denominator=0;
for j=2:length(Time)
  IndexTimeLessJ=Time<Time(j);  %I(yi<yj)
  IndexTimeLargerJ=Time>Time(j); %I(yj<yi)
  IndexLinearPredictLessJ=LinPredictor<LinPredictor(j); % I(betax_j>betax_i)
  IndexLinearPredictLargerJ=LinPredictor>LinPredictor(j); % I(betax_i>betax_j)
  IndexLessThanJ=(1:length(Time))'<j;
  temp=(IndexTimeLessJ & IndexLinearPredictLargerJ & Event);
  if Event(j)
    temp=temp |(IndexTimeLargerJ & IndexLinearPredictLessJ);
  end
  temp=temp & IndexLessThanJ;
  Numerator=Numerator+sum(temp);
  temp=(IndexTimeLessJ & Event);
  if Event(j)
    temp=temp |(IndexTimeLargerJ);
  end
  temp=temp & IndexLessThanJ;
  Denominator=Denominator+sum(temp);
end
I=Numerator/Denominator;
end
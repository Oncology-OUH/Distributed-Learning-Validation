function[Likelihood]=PartialLikelihoodCox(Surv,xval,betaval)
%This code calculated the partial log likelihood for a Cox model. The
%reported number is loglike(null model)-loglike(actual) model. Since the
%likelihood ideally should be larger, so should the log-likelihood. Thus
%the reported number should ideally be a large negative number.
if istable(xval)
  xval=table2array(xval);
end
if istable(betaval)
  betaval=table2array(betaval);
end
if size(betaval,1)==1
  betaval=betaval';
end
if size(betaval,1)~=size(xval,2)
  warning('Size of predictors and betavalues do not match');
end
[~,ix]=sortrows([Surv.Time Surv.Event],[1 2],{'ascend' 'descend'});
Surv=Surv(ix,:);
xval=xval(ix,:);
xbeta=xval*betaval;
theta=exp(xbeta); %This correspojnds to line x in the document
xptheta=diag(theta)*xval;
xpxqtheta=zeros(size(xval,1),length(betaval),length(betaval));
for i=1:size(xval,1)
  temp=xval(i,:);
  xpxqtheta(i,:,:)=theta(i)*(temp')*temp;
end
indexAllEvent=Surv.Event==1;

%Pre-calcluate sum the values at risk for the likelihood function
[~,ia,~]=unique(Surv.Time,'first');
index=Surv.Event(ia)==1;
sumthetaRisk=cumsum(theta,'reverse');
sumthetaRisk=sumthetaRisk(ia);
sumthetaRisk=sumthetaRisk(index);

sumxpthetaRisk=cumsum(xptheta,1,'reverse');
sumxpthetaRisk=sumxpthetaRisk(ia,:);
sumxpthetaRisk=sumxpthetaRisk(index,:);

sumxpxqthetaRisk=cumsum(xpxqtheta,1,'reverse');
sumxpxqthetaRisk=sumxpxqthetaRisk(ia,:,:);
sumxpxqthetaRisk=sumxpxqthetaRisk(index,:,:);


%Pre-calcluate sum the event values for the likelihood function
[uniqueEventTimes,ia,~]=unique(Surv.Time(indexAllEvent),'first');
sumthetaEvents=cumsum(theta(indexAllEvent),'reverse');
sumthetaEvents=sumthetaEvents(ia);
temp=[sumthetaEvents(2:end);0];
sumthetaEvents=sumthetaEvents-temp;

sumxpthetaEvents=cumsum(xptheta(indexAllEvent,:),1,'reverse');
sumxpthetaEvents=sumxpthetaEvents(ia,:);
temp=[sumxpthetaEvents(2:end,:);zeros(1,size(xval,2))];
sumxpthetaEvents=sumxpthetaEvents-temp;

sumxpxqthetaEvents=cumsum(xpxqtheta(indexAllEvent,:,:),1,'reverse');
sumxpxqthetaEvents=sumxpxqthetaEvents(ia,:,:);
temp=[sumxpxqthetaEvents(2:end,:,:);zeros(1,size(xval,2),size(xval,2))];
sumxpxqthetaEvents=sumxpxqthetaEvents-temp;

NrEvents=cumsum(Surv.Event(indexAllEvent),'reverse');
NrEvents=NrEvents(ia);
temp=[NrEvents(2:end);0];
NrEvents=NrEvents-temp;

sumxpbetaEvents=cumsum(xbeta(indexAllEvent),'reverse');
sumxpbetaEvents=sumxpbetaEvents(ia);
temp=[sumxpbetaEvents(2:end);0];
sumxpbetaEvents=sumxpbetaEvents-temp;

sumxpEvents=cumsum(xval(indexAllEvent,:),'reverse');
sumxpEvents=sumxpEvents(ia,:);
temp=[sumxpEvents(2:end,:); zeros(1,size(xval,2))];
sumxpEvents=sumxpEvents-temp;

logLikelihood=0;
gradlogLikelihood=zeros(1,length(betaval));
hesianlogLikelihood=zeros(length(betaval),length(betaval));

%Likelihood.eventTimes=uniqueEventTimes;
%Likelihood.nrEvents=NrEvents;
%Likelihood.sumthetaRisk=sumthetaRisk;
%Likelihood.sumxpthetaRisk=sumxpthetaRisk;
%Likelihood.sumxpxqthetaRisk=sumxpxqthetaRisk;




for j=1:length(uniqueEventTimes)
  logLikelihoodPart=0;
  gradlogLikelihoodPart=zeros(1,length(betaval));
  hessianlogLikelihoodPart=zeros(length(betaval),length(betaval));
  for l=0:(NrEvents(j)-1)
    phijlm=sumthetaRisk(j)-(l/NrEvents(j))*sumthetaEvents(j);
    logLikelihoodPart=logLikelihoodPart+log(phijlm);
    gradphijlm=sumxpthetaRisk(j,:)-(l/NrEvents(j))*sumxpthetaEvents(j,:);
    gradlogLikelihoodPart=gradlogLikelihoodPart+gradphijlm/phijlm;
    doublederivativephijlm=reshape(sumxpxqthetaRisk(j,:,:),size(xval,2),size(xval,2))-(l/NrEvents(j))*reshape(sumxpxqthetaEvents(j,:,:),size(xval,2),size(xval,2));
    hessianlogLikelihoodPart=hessianlogLikelihoodPart+(doublederivativephijlm/phijlm)-(((gradphijlm')*gradphijlm)/(phijlm^2));
  end
  logLikelihood=logLikelihood+sumxpbetaEvents(j)-logLikelihoodPart;
  gradlogLikelihood=gradlogLikelihood+sumxpEvents(j,:)-gradlogLikelihoodPart;
  hesianlogLikelihood=hesianlogLikelihood-hessianlogLikelihoodPart;
end
Likelihood.loglik=logLikelihood;
Likelihood.gradient=gradlogLikelihood';
Likelihood.hessian=hesianlogLikelihood;
end

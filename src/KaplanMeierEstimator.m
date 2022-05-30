function [f,x,flo,fup]=KaplanMeierEstimator(ptdataPIBin,alphaval)
censored=~ptdataPIBin.Event;
[f,x,flo,fup]=ecdf(ptdataPIBin.Time,'Function','survivor','censoring',censored,'Alpha',alphaval);
%The output x includes the minimum value of y as its first two values

if x(1)~=0
  x=[0; x(2:end)];
  f=[1; f(2:end)];
  flo=[1; flo(2:end)];
  fup=[1; fup(2:end)];
else
  x=x(2:end);
  f=f(2:end);
  flo=flo(2:end);
  fup=fup(2:end);
end

if max(x)<max(ptdataPIBin.Time)
  %Need to add censored value to KaplanMeier curve
  x(end+1)=max(ptdataPIBin.Time); 
  f(end+1)=f(end); 
  flo(end+1)=flo(end); 
  fup(end+1)=fup(end); 
end
end
function combinedCalibrationCHarrellKaplanMeierPlots(Data,SavePath)


f1 = figure('Color','w');

for i=1:length(Data)
  subplot(1,length(Data),i);
  yval=Data{i}.CHarrell.MedianBeta.CHarrellBootMedian;
  
  errdown=yval-Data{i}.CHarrell.MedianBeta.CHarrellBootLower;
  errup=Data{i}.CHarrell.MedianBeta.CHarrellBootUpper-yval;
  errorbar(1, yval,errdown, errup,'-ob');
  hold on
  yval=Data{i}.CHarrell.MedianBeta.CV_CHarrellBootMedian;
  errdown=yval-Data{i}.CHarrell.MedianBeta.CV_CHarrellBootLower;
  errup=Data{i}.CHarrell.MedianBeta.CV_CHarrellBootUpper-yval;
  errorbar(2, yval,errdown, errup,'-ob');
  xlim([0.5 2.5]);
  xticks([1 2]);
  xticklabels({'C-Harrell','CV C-Harrell'});
  title(['Center nr: ' num2str(i)]);
end

f2 = figure('Color','w');
for j=1:length(Data)
  
  subplot(1,length(Data),j)
  plot([Data{j}.PI.MedianBeta.CumHist_Levels],[Data{j}.PI.MedianBeta.CumHist],'b');
  hold on
  plot([Data{j}.PI.MedianBeta.CumHist_Levels],[Data{j}.PI.MedianBeta.CumHist_Lower],'b--');
  
  plot([Data{j}.PI.MedianBeta.CumHist_Levels],[Data{j}.PI.MedianBeta.CumHist_Upper],'b--');
  
  ylabel('Cumulative sum');
  
  xlabel('PI');
  title(['Center nr: ' num2str(j)]);
end



f3 = figure('Color','w');
for i=1:length(Data)
  subplot(1,length(Data),i)
  Data{i}.BaselineHazard.MedianBeta=struct2table(Data{i}.BaselineHazard.MedianBeta);
  plot(Data{i}.BaselineHazard.MedianBeta.Time,Data{i}.BaselineHazard.MedianBeta.CumHazard,'r-');
  hold on;
  plot(Data{i}.BaselineHazard.MedianBeta.Time,Data{i}.BaselineHazard.MedianBeta.Lower,'b--');
  plot(Data{i}.BaselineHazard.MedianBeta.Time,Data{i}.BaselineHazard.MedianBeta.Upper,'b--');
  xlabel('Time')
  ylabel('Cumulative Baseline Hazard');
  title(['Center nr: ' num2str(i)]);
end


f4 = figure('Color','w');
for j=1:length(Data)
  for i=1:length(Data{j}.KM_Cox.MedianBeta)
    subplot(1,length(Data),j)
    hold on
    if size(Data{j}.KM_Cox.MedianBeta(i).Values,1)>0
      test=Data{j}.KM_Cox.MedianBeta(i).Values;
      test=struct2cell(test);
      index=cellfun(@(x) isempty(x),test);
      test(index)={NaN};
      Data{j}.KM_Cox.MedianBeta(i).Values=cell2table(test','VariableNames',fieldnames(Data{j}.KM_Cox.MedianBeta(i).Values));
      %Data{j}.KM_Cox.MedianBeta(i).Values=struct2table(Data{j}.KM_Cox.MedianBeta(i).Values);
      plot(Data{j}.KM_Cox.MedianBeta(i).Values.Time,Data{j}.KM_Cox.MedianBeta(i).Values.KM,'r');
      plot(Data{j}.KM_Cox.MedianBeta(i).Values.Time,Data{j}.KM_Cox.MedianBeta(i).Values.KM_Lower,'b--');
      plot(Data{j}.KM_Cox.MedianBeta(i).Values.Time,Data{j}.KM_Cox.MedianBeta(i).Values.KM_Upper,'b--');
      plot(Data{j}.KM_Cox.MedianBeta(i).Values.Time,Data{j}.KM_Cox.MedianBeta(i).Values.Cox,'g-');
    end
    xlabel('Time');
    ylabel('Survival - KM(red) - Cox(green)');
  end
  title(['Center nr: ' num2str(j)]);
end


f5 = figure('Color','w');
for j=1:length(Data)
  for i=1:length(Data{j}.CalPlot.MedianBeta)
    subplot(length(Data{j}.CalPlot.MedianBeta),length(Data),(i-1)*length(Data)+j);
    hold on
    xval=[Data{j}.CalPlot.MedianBeta(i).Plot.Predicted];
    yval=[Data{j}.CalPlot.MedianBeta(i).Plot.KM_Observed];
    elow=yval - [Data{j}.CalPlot.MedianBeta(i).Plot.Lower];
    eup=[Data{j}.CalPlot.MedianBeta(i).Plot.Upper]-yval;
    errorbar(xval,yval,elow,eup,'o');
    plot([0,1],[0,1],'b-')
    xlabel(['Predicted surival at time: ',num2str(Data{j}.CalPlot.MedianBeta(i).time)]);
    ylabel('Observed surival');
    if i==1
      title(['Center nr: ' num2str(j)]);
    end
  end
end
tempdate=datestr(datetime);
tempdate=strrep(tempdate,':','_');
tempdate=strrep(tempdate,' ','_');
temp=[fullfile(SavePath,'CHarrellPlot'),'_',tempdate];
savefig(f1,temp);
temp=[fullfile(SavePath,'CumHistPI'),'_',tempdate];
savefig(f2,temp)
temp=[fullfile(SavePath,'BaselineHazardPlot'),'_',tempdate];
savefig(f3,temp);
temp=[fullfile(SavePath,'KMPlot'),'_',tempdate];
savefig(f4,temp);
temp=[fullfile(SavePath,'CalibrationPlot'),'_',tempdate];
savefig(f5,temp);

end
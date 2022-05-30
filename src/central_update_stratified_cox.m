function [model]=central_update_stratified_cox(model)

tollerance=model.ConvergTollerance;
NBoot=size(model.beta(1).betaval,1);
NModels=length(model.beta);
NrCenter=length(model.client);
Converted=false(NModels,1);
for k=1:NModels
  hessian=cell(NBoot,1);
  hessian(:)={zeros(size(model.client(1).I(k).Result(1).hessian))};
  gradient=cell(NBoot,1);
  gradient(:)={zeros(size(model.client(1).I(k).Result(1).gradient))};

  for j=1:NrCenter
    gradient=cellfun(@(x,y) x+y,gradient,{model.client(j).I(k).Result.gradient}','UniformOutput',false);
    hessian=cellfun(@(x,y) x+y,hessian,{model.client(j).I(k).Result.hessian}','UniformOutput',false);
  end
  deltabeta=cellfun(@(x,y) (x\(-y))',hessian,gradient,'UniformOutput',false);
  belowtol_index=cellfun(@(x) all(abs(x)<tollerance),deltabeta);
  deltabeta(belowtol_index)={zeros(size(deltabeta{1}))};
  tempvarnames=model.client(1).I(k).VarNames;
  betacurrent=model.beta(k).betaval(:,tempvarnames);  
  model.beta(k).betaval=array2table(table2array(betacurrent)+cell2mat(deltabeta),'VariableName',tempvarnames);
  
  if all(belowtol_index)
    Converted(k)=true;
  end
end
if all(Converted)
    model.convergence=true;
end

% 
% Nboot = length(model.randseed);
% NModel = size(model.varselect,1);
% gradient_all=zeros(size(model.client(1).gradient));
% hessian_all=zeros(size(model.client(1).hessian));
% beta_update=zeros(Nboot*NModel,size(model.varselect,2));
% single_beta_update = zeros(1,size(model.varselect,2));
% for j=1:length(model.client)
%   gradient_all = gradient_all+model.client(j).gradient;
%   hessian_all = hessian_all+model.client(j).hessian;
% end
% 
% for j=1:NModel
%   varind = logical(model.varselect(j,:));
%   varlen = sum(model.varselect(j,:));
%   for i=1:Nboot
%     deltabeta = (hessian_all((j-1)*Nboot+i,1:varlen^2) \ (-gradient_all((j-1)*Nboot+i,1:varlen)));
%     beta_update((j-1)*Nboot+i,varind) = deltabeta(1,:);
%     single_beta_update(varind) = deltabeta(1,:);
%     model.beta(j).betaval(i,:) = model.beta(j).betaval(i,:)+single_beta_update;
%   end
% end
% 
% if all(beta_update(:)' < model.toler)
%   model.convergence=true;
% end

end
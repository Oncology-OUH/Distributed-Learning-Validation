function LocalLarynxValidation(MPI_config,MPI_name)
%Read the local config file
config = readConfigFile(MPI_config);

%Except for those variables contained in DoNotMoveToModel below, transfer
%all local config variables to the model variable. All string parameters
%that can be converted to a number are converted. This change is made to
%facilitate easy implementation of additional parameters to the local
%config file.
varnames=fieldnames(config);
DoNotMoveToModel={'dataLocation','MPIClientPort'};
for i=1:length(varnames)
  if sum(strcmp(varnames{i},DoNotMoveToModel))==0
    model.(varnames{i})=config.(varnames{i});
    temp=str2double(model.(varnames{i}));
    if ~isnan(temp)
      model.(varnames{i})=temp;
    end
  end
end

%The next line initializes the local Message Passing Interface (MPI), wait
%until the central server report that it is ready, and finally add the
%local node to the optimization network. If successful mpi.startedMPI is
%returned as true.
mpi = initialise_client_comm(config,MPI_name);


if mpi.startedMPI
    prev_itr=-1;
    model.state='';
    while (~strcmp(model.state,'terminate_function')) && (mpi.timeouts<5)
        
        model.alg_client_index=1;
        %model.client_name = config.clientName;
        %model.NrPtPerBin=str2double(config.NrPtPerBin);
        %model.localseed=str2double(config.seed);
        
        [itr, model, ~, central_message, mpi] = Client_Receive_from_MPI(mpi,prev_itr,model);
        [model, ~] = readJSONMessage('local',model,central_message,model.alg_client_index);
        
        [data, ~,model]= loadLarynxData(config, model);
        % Partition data for validation procedures
        
        %[model, data] = initialize_validations(model, data, model.validation);
                
        while (~strcmp(model.state,'terminate_algorithm')) && (mpi.timeouts<5)
            [var,model,mpi,itr]=local_stratified_cox(mpi,itr,central_message,model,data);
            [itr,model,central_message, mpi]=sendMessage(mpi,model,itr,mpi.ClientIndex,var);
        end
        disp('Terminating algorithm session and re-initialising.')
        [itr,model,~,mpi]=sendMessage(mpi,model,itr,mpi.ClientIndex,[]);
    end
    disp('Terminating function.')
    sendReplyToMaster(mpi.svc, itr, 'null');     
end

end


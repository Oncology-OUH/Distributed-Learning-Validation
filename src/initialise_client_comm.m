function mpi = initialise_client_comm(config,MPI_name)

%Define the location of the local Message Passing Interface (MPI). It is
%the local MPI that communicate with the central communication server
mpi=struct;
if ~isfield(config, 'MPIClientPort')
  config.MPIClientPort='80';
end
mpi.svc=MPIService(['http://localhost:' config.MPIClientPort '/' MPI_name '/MPIClient'],...
  ['http://localhost:' config.MPIClientPort '/' MPI_name '/MPIClient?wsdl']);

disp(['Client starting from ' config.client_name])

%Initialize the local MPI (Message Passing Interface) client that is used
%to communicate with the central communication server (MPI server). If the
%initialization is okay true is returned from the fuction

mpi.startedMPI = initialiseClientMPI(mpi.svc);

mpi.timeouts=0;
if strcmp(mpi.startedMPI,'false')
  mpi.startedMPI=false;
  mpi.comm_message = 'MPI failed to start'; disp(mpi.comm_message)
else
  %The local MPI is now initialized. Now, wait for a response that the
  %central server is ready
  disp('Initialisation of MPI successful. Waiting for central algorithm to start');
  while strcmp(isMasterOn(mpi.svc),'false');  pause(0.25);  end
  disp('Central algorithm is ON');
  %Inform the MPI system that the local client is ready, and report whether the addition of the local client was successful
  mpi.ClientIndex=str2double(addClient(mpi.svc));
  if mpi.ClientIndex==-1
    mpi.comm_message = ['Cannot add Client: ' config.client_name]; disp(mpi.comm_message)
    mpi.startedMPI = false;
  else
    mpi.comm_message = strcat(['Client ' config.client_name ' added with index: '],num2str(mpi.ClientIndex)); disp(mpi.comm_message)
    mpi.startedMPI = true;
  end
  
end

end

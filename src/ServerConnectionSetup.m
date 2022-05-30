function MPI = ServerConnectionSetup(MPI_name,participating_clients,initialise,PortNumber)
%Have included the option to send the port number to the function. If the variable is not provided, the port number will default to port 80.
if ~exist('PortNumber','var')
  PortNumber='80';
end
%%% Process input arguments
number_of_participating_clients=length(strsplit(participating_clients,','));
number_of_clients=0;
disp(['Number of participating clients: ' num2str(number_of_participating_clients)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Connect to distribution server
svc=MPIService(['http://localhost:', PortNumber ,'/' MPI_name '/MPIClient'],['http://localhost:',PortNumber,'/' MPI_name '/MPIClient?wsdl']);
startedMPI = initialiseCentralMPI(svc);
if strcmp(startedMPI,'false')
    disp('MPI failed to start')
else
    disp('Initialisation of MPI successful')
    if initialise
        initClientInfo(svc);
    end
    
    if str2double(getNumberofCurrentClients(svc))~=number_of_participating_clients
        setAllowedClientsList(svc,participating_clients);
        setMasterOn(svc);
        
        %%% Wait for each client to also connect to the server
        while str2double(getNumberofCurrentClients(svc))<number_of_participating_clients
%             [str2double(getNumberofCurrentClients(svc)) number_of_participating_clients]
            
            pause(0.5);
        end
    end
    number_of_clients= str2double(getNumberofCurrentClients(svc));
    
    disp('All clients have connected to the server')
    
end
MPI.svc=svc;
MPI.num_clients = number_of_clients;

end
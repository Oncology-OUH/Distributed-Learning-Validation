function [model, itr] = Broadcast_Receive_from_MPI(svc,master_message,num_clients,itr,model)
%This function sends a request and parameter data to all the local nodes.
%The actual message is sent to the local nodes using the command
%publishMessageToClients. After having received the request locally, the
%local nodes process the request based on the supplied parameters and the
%local data (only available at the local node) and send a message back to
%the MPI system. Whether all local nodes have replied to the request is
%checked in a loop every 0.2 seconds. After having detected that all nodes
%have replied, all the replies are collected from the MPI system and
%returned
itr=itr+1;
publishMessageToClients(svc,master_message,itr); % Send the model update to all clients in network

retry = 0; success=false;
while retry < 5 && ~success
    try
        %%%Wait until all clients have received model
        while strcmp(didAllClientsReply(svc,itr),'false')
            pause(0.2);
        end
        success=true;
    catch err
        retry=retry+1;
        disp(['Encountered an exception relating to function call "didAllClientsReply()". Attempt: ' num2str(retry)])
        disp(err.message);
        for st = 1:size(err.stack,1)
            disp(err.stack(st).file);disp(err.stack(st).name);disp(err.stack(st).line);
        end  
    end
end
if ~success
    error(['Error after ' num2str(retry) ' attempts recovering from exceptions thrown from "didAllClientsReply()".'])
end


for i=1:num_clients   % Collect the updated model from each client
    
    retry = 0; success=false;
    while retry < 5 && ~success
        try
            client_msg=getClientReply(svc,i-1);
            success=true;
        catch err
            retry=retry+1;
            disp(['Encountered an exception relating to function call "getClientReply()" for client: ' num2str(i) '. Attempt: ' num2str(retry)])
            disp(err.message);
            for st = 1:size(err.stack,1)
                disp(err.stack(st).file);disp(err.stack(st).name);disp(err.stack(st).line);
            end
        end
    end
    if ~success
        error(['Error after ' num2str(retry) ' attempts recovering from exceptions thrown from "getClientReply()" for client: ' num2str(i)])
    end
    
    [model, received_parameters] = readJSONMessage('central',model,client_msg,i);
end


end
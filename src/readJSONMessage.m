
function [model, received_parameters] = readJSONMessage(mode,model,message,client_index)

    if strcmp(mode,'central')
        temp = jsondecode(message);
        if ~isempty(temp)
            msg_fields = fieldnames(temp);
            for i=1:length(msg_fields)
                model.client(client_index).(msg_fields{i}) = temp.(msg_fields{i});
            end
            received_parameters = fieldnames(model.client(client_index));
        else
            received_parameters=[];
        end
        
    elseif strcmp(mode,'local')
    
        model.central = jsondecode(message);
        received_parameters = fieldnames(model.central);
        for j=1:length(received_parameters)
            model.(received_parameters{j}) = model.central.(received_parameters{j});
        end
    end

end
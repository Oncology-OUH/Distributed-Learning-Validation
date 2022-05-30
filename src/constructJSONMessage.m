function message = constructJSONMessage(model,parameters)
    packet=[];
    for i=1:length(parameters)
        packet.(parameters{i}) = model.(parameters{i});
    end
    message = jsonencode(packet);

end

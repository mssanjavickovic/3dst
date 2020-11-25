function [All_new_data, new_channel_name,opt_gate_context] = Merge_ST_Clusters(CSV_name,All_new_data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Retrive data from histoCAT
handles = gethand;
selected_gates = get(handles.list_samples,'Value');
gates = retr('gates');
sessionData = retr('sessionData');
gate_context = retr('gateContext');

gates_names_aligned = zeros(size(gates,1),2);
gates_names_split = [];
gates_names_cell = {};

for i=1:size(gates,1)
    % Extract information from gates
    gates_names_split{i} = strsplit(gates{i,1},'_');
    gates_names_aligned(i,1) = str2num(gates_names_split{i}{1,4});
    gates_names_aligned(i,2) = str2num(gates_names_split{i}{1,5});
    gates_names_cell{i} = strcat(gates_names_split{i}{1,4},'_',...
        gates_names_split{i}{1,5});
end

% Load custers
Cluster_load_raw = readtable(CSV_name);
Cluster_aligned = table2array(Cluster_load_raw(:,2:3));

Cluster_names_cell = {};
Find_cluster_location = [];

for k=1:size(Cluster_aligned,1)
    Cluster_names_cell{k} = strcat(num2str(Cluster_aligned{k,1}),'_',...
        num2str(Cluster_aligned{k,2}));
    Find_cluster_location(k,1) = find(strcmp...
        (gates_names_cell,Cluster_names_cell{k}));

    
    %create new session data
    opt_gates = Find_cluster_location(k,1);
    opt_gate_context = gates{opt_gates,2};
        
    % New channel names
    new_channel_name = Cluster_load_raw.Properties.VariableNames(4:16);
    
    % Extract info from CSV file
    new_data_oneline = table2array(Cluster_load_raw(k,4:16));
    new_data = repmat(new_data_oneline,size(gates{opt_gates,2},2),1);
    
    All_new_data(opt_gate_context,1:size(new_data,2)) = new_data;
    
end

end


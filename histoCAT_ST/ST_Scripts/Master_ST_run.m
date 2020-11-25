% Retrive data from histoCAT
handles = gethand;
selected_gates = get(handles.list_samples,'Value');
gates = retr('gates');
sessionData = retr('sessionData');
gate_context = retr('gateContext');
opt_gates = 1:size(gates,1);
opt_gate_context = 1:size(sessionData,1);

%CSV file
RA_which = 'RA2';
Sample_which = '7';

% Load Cell Types
CSV_name = strcat('/Users/denis/Desktop/ST/CellTypes/Cell types/RA2/','scRNAseq.',Sample_which,'.',RA_which,'.csv');
Cell_type_raw = readtable(CSV_name);
%Cluster_aligned = table2array(Cluster_load_raw(:,2:3));
Amount_markers = length(4:16);

% Create zero array for all new data
All_new_data = zeros(size(sessionData,1),Amount_markers);

%Cell Types
[All_new_data,new_channel_name] = Merge_ST_Clusters(strcat('/Users/denis/Desktop/ST/CellTypes/Cell types/RA2/','scRNAseq.',Sample_which,'.',RA_which,'.csv'),All_new_data);

[All_new_data,~] = Merge_ST_Clusters(strcat('/Users/denis/Desktop/ST/','Cluster1.',Sample_which,'.',RA_which,'.csv'),All_new_data);
[All_new_data,~] = Merge_ST_Clusters(strcat('/Users/denis/Desktop/ST/','Cluster2.',Sample_which,'.',RA_which,'.csv'),All_new_data);
[All_new_data,~] = Merge_ST_Clusters(strcat('/Users/denis/Desktop/ST/','Cluster3.',Sample_which,'.',RA_which,'.csv'),All_new_data);
[All_new_data,new_channel_name] = Merge_ST_Clusters(strcat('/Users/denis/Desktop/ST/','Cluster4.',Sample_which,'.',RA_which,'.csv'),All_new_data);

%[All_new_data,new_channel_name] = Merge_ST_Clusters_batch(strcat('/Users/denis/Desktop/ST/scRNAseq.2.RA1.csv'),All_new_data,Cluster_load_raw,Cluster_aligned);

% Write session Data
 [sessionData,gates] =  addChannels(new_channel_name,...
     All_new_data, opt_gate_context,opt_gates, gates, sessionData);

[pheno_clusters] = parse_Phenographclusters(sessionData,gates);
custom_gatesfolder = retr('custom_gatesfolder');





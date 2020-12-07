%Create tSNE and PhenoGraph from ST RA data
RA2_table = table();
Names = {'RA2_1','RA2_2','RA2_3','RA2_4','RA2_5','RA2_6','RA2_7'};

for i=1:size(Names,2)
    % Create temporary table for each RA section
    %RA2_table_temp = readtable(strcat(Names{1,i},'.csv'));
    RA2_table_temp = readtable(strcat(Names{1,i},'_cellTypes.csv'));
    % Find neighbors variables (which make table unmergeable)
    hasMatch = cellfun('isempty',...
        regexp(RA2_table_temp.Properties.VariableNames, 'neighbour', 'once'));
    
    hasMatch_PhenotypicCluster = cellfun('isempty',...
        regexp(RA2_table_temp.Properties.VariableNames, 'PhenotypicCluster', 'once'));
    
    hasMatch(1,find(~hasMatch_PhenotypicCluster)) = 0;
    
    
%     % Find cellTypes
%     hasMatch = cellfun('isempty',...
%         regexp(RA2_table_temp.Properties.VariableNames, 'neighbour', 'once'));
%     
    % Remove neighbors to enable merging
    RA2_table_temp_short = RA2_table_temp...
        (:, RA2_table_temp.Properties.VariableNames(hasMatch));
    if i == 1
        % Create for first round
        RA2_table =RA2_table_temp_short;
        
    else
        % Merge
        RA2_table = [RA2_table;RA2_table_temp_short];
        RA2_table_temp_short = [];
    end
end

% Run_tSNE
% Select which channels to use for tSNE
data = table2array(RA2_table(:,([3:7 17])));
% normalize
percentile = 99;
data = mynormalize(data, percentile);
% Run tSNE
try
    map = fast_tsne(data, 110);
catch
    msgbox(...
        ['tSNE Failed: Common causes are \n' ...
        'a) illegal histoCAT installation path - spaces in path.\n' ...
        'b) illegal histoCAT installation path - no writing persmissions in folder.\n' ...
        'c) perplexity too high caused by insufficient number of points.'],...
        'Error','error');
    return;
end

%save('tSNE_all_RA2_new.mat')
load('tSNE_all_RA2_new.mat')

% Visualize
RA2_table.tSNE1 = map(:,1);
RA2_table.tSNE2 = map(:,2);

% Remove large cells
remove_large_area = RA2_table.Area > 3000;
remove_large_area_index = find(remove_large_area);
RA2_table.Area(remove_large_area_index) = 3000;

figure()
scatter(RA2_table.tSNE1,RA2_table.tSNE2,[],...
    RA2_table.Percent_Touching,...
    'LineWidth',0.5);

% Run PhenoGraph
[labels,~,~] = phenograph(data, 150,'random_seed','No');

figure()
scatter(RA2_table.tSNE1,RA2_table.tSNE2,[],...
    labels,...
    'LineWidth',0.5);


% Add more cell types
RA2_table.macrophage_clusters = RA2_table_celltypes.macrophage_clusters;

RA2_table.fibro2a_thy1 = RA2_table_celltypes.fibro2a_thy1;
RA2_table.fibro1_cd55 = RA2_table_celltypes.fibro1_cd55;
RA2_table.fibro2b_thy1 = RA2_table_celltypes.fibro2b_thy1;

figure()
scatter(RA2_table.tSNE1,RA2_table.tSNE2,[],...
    RA2_table.macrophage_clusters,...
    'LineWidth',0.5);
title('macrophage_clusters')

figure()
scatter(RA2_table.tSNE1,RA2_table.tSNE2,[],...
    RA2_table.fibro1_cd55,...
    'LineWidth',0.5);
title('fibro1_cd55')

figure()
scatter(RA2_table.tSNE1,RA2_table.tSNE2,[],...
    RA2_table.fibro2a_thy1,...
    'LineWidth',0.5);
title('fibro2a_thy1')

figure()
scatter(RA2_table.tSNE1,RA2_table.tSNE2,[],...
    RA2_table.fibro2b_thy1,...
    'LineWidth',0.5);
title('fibro2b_thy1')


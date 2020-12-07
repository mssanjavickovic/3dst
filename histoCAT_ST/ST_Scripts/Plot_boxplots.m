% Run with histoCAT open
load('RA1_1_noCluster.mat');
[Average_table_RA1_1] = ST_Calculate_Density('RA1_1',gates);
load('RA1_2_session_withNoCluster.mat');
[Average_table_RA1_2] = ST_Calculate_Density('RA1_2',gates);
load('sessionData_RA1_3.mat');
[Average_table_RA1_3] = ST_Calculate_Density('RA1_3',gates);
load('RA1_4_sessionData.mat');
[Average_table_RA1_4] = ST_Calculate_Density('RA1_4',gates);

% Create table for boxplot all together
Average_table_plot(:,1) = [Average_table_RA1_1.meanPercentTouching;...
    Average_table_RA1_2.meanPercentTouching;...
    Average_table_RA1_3.meanPercentTouching;...
    Average_table_RA1_4.meanPercentTouching];
% % Create a second column for cluster
% Average_table_plot(:,2) = [Average_table_RA1_1.cluster;...
%     Average_table_RA1_2.cluster;...
%     Average_table_RA1_3.cluster;...
%     Average_table_RA1_4.cluster];
% % Where is cluster one?
% Cluster_one_logic = (Average_table_plot(:,2) == 1);
% Cluster_one_index = find(Cluster_one_logic);
% % Select only cluster 1 density
% Density_cluster = Average_table_plot(Cluster_one_logic,1)>0.7;
% Density_percent = (size(find(Density_cluster),1))/(size(Density_cluster,1));



%boxplot
g = [ones(size(Average_table_RA1_1.meanPercentTouching));...
    2*ones(size(Average_table_RA1_2.meanPercentTouching));...
    3*ones(size(Average_table_RA1_3.meanPercentTouching));...
    4*ones(size(Average_table_RA1_4.meanPercentTouching))];

figure
boxplot(Average_table_plot,g)
Names = {'RA1_1','RA1_2','RA1_3','RA1_4'};
set(gca,'XTickLabel',Names);
ylabel('Percent Touching')

hold on
rng(2)
gscatter(ones(size(Average_table_RA1_1.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA1_1.meanPercentTouching))-0.5)/2)...
    ,Average_table_RA1_1.meanPercentTouching,Average_table_RA1_1.inf,'br','xo');
hold on
rng(2)
gscatter(2*ones(size(Average_table_RA1_2.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA1_2.meanPercentTouching))-0.5)/4)...
    ,Average_table_RA1_2.meanPercentTouching,Average_table_RA1_2.inf,'br','xo');
hold on
rng(2)
gscatter(3*ones(size(Average_table_RA1_3.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA1_3.meanPercentTouching))-0.5)/6)...
    ,Average_table_RA1_3.meanPercentTouching,Average_table_RA1_3.inf,'br','xo');
hold on
rng(2)
gscatter(4*ones(size(Average_table_RA1_4.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA1_4.meanPercentTouching))-0.5)/8)...
    ,Average_table_RA1_4.meanPercentTouching,Average_table_RA1_4.inf,'br','xo');

figure
boxplot(Average_table_plot,g)
Names = {'RA1_1','RA1_2','RA1_3','RA1_4'};
set(gca,'XTickLabel',Names);
ylabel('Percent Touching')

hold on
rng(2)
gscatter(ones(size(Average_table_RA1_1.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA1_1.meanPercentTouching))-0.5)/2)...
    ,Average_table_RA1_1.meanPercentTouching,Average_table_RA1_1.cluster);
hold on
rng(2)
gscatter(2*ones(size(Average_table_RA1_2.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA1_2.meanPercentTouching))-0.5)/4)...
    ,Average_table_RA1_2.meanPercentTouching,Average_table_RA1_2.cluster);
hold on
rng(2)
gscatter(3*ones(size(Average_table_RA1_3.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA1_3.meanPercentTouching))-0.5)/6)...
    ,Average_table_RA1_3.meanPercentTouching,Average_table_RA1_3.cluster);
hold on
rng(2)
gscatter(4*ones(size(Average_table_RA1_4.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA1_4.meanPercentTouching))-0.5)/8)...
    ,Average_table_RA1_4.meanPercentTouching,Average_table_RA1_4.cluster);


% Run for RA2
% Run with histoCAT open
load('RA2_1.mat');
[Average_table_RA2_1] = ST_Calculate_Density('RA2_1',gates);
clear gates

load('RA2_2.mat');
[Average_table_RA2_2] = ST_Calculate_Density('RA2_2',gates);
clear gates

load('RA2_3.mat');
[Average_table_RA2_3] = ST_Calculate_Density('RA2_3',gates);
clear gates

load('RA2_4.mat');
[Average_table_RA2_4] = ST_Calculate_Density('RA2_4',gates);
clear gates

load('RA2_5.mat');
[Average_table_RA2_5] = ST_Calculate_Density('RA2_5',gates);
clear gates

load('RA2_6.mat');
[Average_table_RA2_6] = ST_Calculate_Density('RA2_6',gates);
clear gates

load('RA2_7.mat');
[Average_table_RA2_7] = ST_Calculate_Density('RA2_7',gates);
clear gates

% load('RA2_8.mat');
% [Average_table_RA2_8] = ST_Calculate_Density('RA2_8',gates);
% clear gates
% 
% load('RA2_9.mat');
% [Average_table_RA2_9] = ST_Calculate_Density('RA2_9',gates);
% clear gates

% Create table for boxplot all together
Average_table_plot_RA2(:,1) = [Average_table_RA2_1.meanPercentTouching;...
    Average_table_RA2_2.meanPercentTouching;...
    Average_table_RA2_3.meanPercentTouching;...
    Average_table_RA2_4.meanPercentTouching;...
    Average_table_RA2_5.meanPercentTouching;...
    Average_table_RA2_6.meanPercentTouching;...
    Average_table_RA2_7.meanPercentTouching];

% Create a second column for cluster
Average_table_plot_RA2(:,2) = [Average_table_RA2_1.cluster;...
    Average_table_RA2_2.cluster;...
    Average_table_RA2_3.cluster;...
    Average_table_RA2_4.cluster;...
    Average_table_RA2_5.cluster;...
    Average_table_RA2_6.cluster;...
    Average_table_RA2_7.cluster;];

% Where is cluster one?
Cluster_one_logic = (Average_table_plot_RA2(:,2) == 1);
Cluster_one_index = find(Cluster_one_logic);
% Select only cluster 1 density
Density_cluster = Average_table_plot_RA2(Cluster_one_logic,1)>0.7;
Density_percent = (size(find(Density_cluster),1))/(size(Density_cluster,1));


%boxplot
g = [ones(size(Average_table_RA2_1.meanPercentTouching));...
    2*ones(size(Average_table_RA2_2.meanPercentTouching));...
    3*ones(size(Average_table_RA2_3.meanPercentTouching));...
    4*ones(size(Average_table_RA2_4.meanPercentTouching));...
    5*ones(size(Average_table_RA2_5.meanPercentTouching));...
    6*ones(size(Average_table_RA2_6.meanPercentTouching));...
    7*ones(size(Average_table_RA2_7.meanPercentTouching));];

figure
boxplot(Average_table_plot_RA2,g)
Names = {'RA2_1','RA2_2','RA2_3','RA2_4','RA2_5','RA2_6','RA2_7'};
set(gca,'XTickLabel',Names);
ylabel('Percent Touching')

hold on
rng(2)
gscatter(ones(size(Average_table_RA2_1.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA2_1.meanPercentTouching))-0.5)/2)...
    ,Average_table_RA2_1.meanPercentTouching,Average_table_RA2_1.inf,'br','xo');
hold on
rng(2)
gscatter(2*ones(size(Average_table_RA2_2.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA2_2.meanPercentTouching))-0.5)/4)...
    ,Average_table_RA2_2.meanPercentTouching,Average_table_RA2_2.inf,'br','xo');
hold on
rng(2)
gscatter(3*ones(size(Average_table_RA2_3.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA2_3.meanPercentTouching))-0.5)/6)...
    ,Average_table_RA2_3.meanPercentTouching,Average_table_RA2_3.inf,'br','xo');
hold on
rng(2)
gscatter(4*ones(size(Average_table_RA2_4.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA2_4.meanPercentTouching))-0.5)/8)...
    ,Average_table_RA2_4.meanPercentTouching,Average_table_RA2_4.inf,'br','xo');
hold on
rng(2)
gscatter(5*ones(size(Average_table_RA2_5.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA2_5.meanPercentTouching))-0.5)/10)...
    ,Average_table_RA2_5.meanPercentTouching,Average_table_RA2_5.inf,'br','xo');
hold on
rng(2)
gscatter(6*ones(size(Average_table_RA2_6.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA2_6.meanPercentTouching))-0.5)/12)...
    ,Average_table_RA2_6.meanPercentTouching,Average_table_RA2_6.inf,'br','xo');
hold on
rng(2)
gscatter(7*ones(size(Average_table_RA2_7.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA2_7.meanPercentTouching))-0.5)/14)...
    ,Average_table_RA2_7.meanPercentTouching,Average_table_RA2_7.inf,'br','xo');

figure
boxplot(Average_table_plot_RA2,g)
Names = {'RA2_1','RA2_2','RA2_3','RA2_4','RA2_5','RA2_6','RA2_7'};
set(gca,'XTickLabel',Names);
ylabel('Percent Touching')

hold on
rng(2)
gscatter(ones(size(Average_table_RA2_1.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA2_1.meanPercentTouching))-0.5)/2)...
    ,Average_table_RA2_1.meanPercentTouching,Average_table_RA2_1.cluster);
hold on
rng(2)
gscatter(2*ones(size(Average_table_RA2_2.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA2_2.meanPercentTouching))-0.5)/4)...
    ,Average_table_RA2_2.meanPercentTouching,Average_table_RA2_2.cluster);
hold on
rng(2)
gscatter(3*ones(size(Average_table_RA2_3.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA2_3.meanPercentTouching))-0.5)/6)...
    ,Average_table_RA2_3.meanPercentTouching,Average_table_RA2_3.cluster);
hold on
rng(2)
gscatter(4*ones(size(Average_table_RA2_4.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA2_4.meanPercentTouching))-0.5)/8)...
    ,Average_table_RA2_4.meanPercentTouching,Average_table_RA2_4.cluster);
hold on
rng(2)
gscatter(5*ones(size(Average_table_RA2_5.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA2_5.meanPercentTouching))-0.5)/10)...
    ,Average_table_RA2_5.meanPercentTouching,Average_table_RA2_5.cluster);
hold on
rng(2)
gscatter(6*ones(size(Average_table_RA2_6.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA2_6.meanPercentTouching))-0.5)/12)...
    ,Average_table_RA2_6.meanPercentTouching,Average_table_RA2_6.cluster);
hold on
rng(2)
gscatter(7*ones(size(Average_table_RA2_7.meanPercentTouching))...
    .*(1+(rand(size(Average_table_RA2_7.meanPercentTouching))-0.5)/14)...
    ,Average_table_RA2_7.meanPercentTouching,Average_table_RA2_7.cluster);
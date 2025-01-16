clc; clear all; close all;

%% Loading Data
path = '/mfip/mfip1/arielle/PhDProject2/';
%parcellation = strsplit(input("What parcellation scheme and subject database do you want to do? (i.e.: schaefer 100 MICs HC) \n" , "s"));
  %  parcelName = 'schaefer'; %parcellation{1};
  %  parcelNum  = str2double(parcellation{2});
    groupName = 'MICs'; %parcellation{3};
  %  subjType = parcellation{4};

data = readtable(fullfile(path, 'results', groupName, 'clusteringCoefficient', 'averageClusteringCoefficients.csv'));

% Converting Parcellations to Categorical values
data.Parcellation = categorical(data.Parcellation);
colors = [[251/255 180/255 174/255]; [254/255 217/255 166/255]; [255/255 238/255 140/255]; [204/255 235/255 197/255]; [179/255 205/255 227/255]; [253/255 218/255 236/255]; ...
    [222/255 203/255 228/255]; [205/255 194/255 170/255]; [227/255 227/255 227/255]]; 

%% Plotting HC
% data but just HCs
data_HC = data(strcmp(data.Group, 'HC'), :);

figure; 
hold on; 
groups_HC = categories(data_HC.Parcellation);

for i = 1:numel(groups_HC)
    % Select the data for the current group
    groupData = data_HC(data_HC.Parcellation == groups_HC{i}, :);
    groupData = sortrows(groupData,'Threshold','ascend');
    % Plot Threshold vs. averageclusteringCoefficient for the current group
    p = plot(groupData.Threshold, groupData.GlobalClusteringCoefficient,'o:', ...
        'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i, :), ...
        'DisplayName', groups_HC{i}, 'Color', colors(i, :), 'LineWidth', 4 );
end

% Add labels and legend
xlabel('Threshold');
ylabel('Average Clustering Coefficient');
ylim([min(data.GlobalClusteringCoefficient) max(data.GlobalClusteringCoefficient)])
legend('show');
title('Threshold vs. Average Clustering Coefficient HC');

% Optional: Adjust axis limits for better visualization
axis tight;
grid on;

% Release hold
hold off;
figure; 


%% Plotting PX
data_PX = data(strcmp(data.Group, 'PX'), :);
groups_PX = categories(data_PX.Parcellation);

hold on;
for i = 1:numel(groups_PX)
    % Select the data for the current group
    groupData = data_PX(data_PX.Parcellation == groups_PX{i}, :);
    groupData = sortrows(groupData,'Threshold','ascend');
    % Plot Threshold vs. averageclusteringCoefficient for the current group
    p = plot(groupData.Threshold, groupData.GlobalClusteringCoefficient,'o:', ...
        'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i, :), ...
        'DisplayName', groups_HC{i}, 'Color', colors(i, :), 'LineWidth', 4 );
end

% Add labels and legend
xlabel('Threshold');
ylabel('Average Clustering Coefficient');
ylim([min(data.GlobalClusteringCoefficient) max(data.GlobalClusteringCoefficient)])
legend('show');
title('Threshold vs. Average Clustering Coefficient PX');

% Optional: Adjust axis limits for better visualization
axis tight;
grid on;

% Release hold
hold off;

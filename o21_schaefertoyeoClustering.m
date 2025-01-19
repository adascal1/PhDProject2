clc; clear all; close all; 

% Load the LUT file
filename = '/mfip/mfip1/arielle/MICs_dataset/MNI/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI/fsleyes_lut/Schaefer2018_900Parcels_7Networks_order.lut';
data = readtable(filename, 'FileType', 'text', 'Delimiter', ' ', 'ReadVariableNames', false);

% Extract relevant columns
indices = data.Var1;         % First column
labels = data.Var5;          % Fifth column (contains the full string)

% Initialize a map for grouping
groupedValues = containers.Map;

for i = 1:length(labels)
    % Split the label by underscores
    parts = split(labels{i}, '_');
    
    % Extract the third part
    if length(parts) >= 3
        key = parts{3};
    else
        key = 'Unknown'; % Fallback for malformed strings
    end
    
    % Add the index to the appropriate group
    if isKey(groupedValues, key)
        groupedValues(key) = [groupedValues(key); indices(i)];
    else
        groupedValues(key) = indices(i);
    end
end

% Display the grouped results
disp('Grouped Values:');
keys = groupedValues.keys;
for i = 1:length(keys)
    fprintf('Group: %s\n', keys{i});
    disp(groupedValues(keys{i}));
end

% Save the grouped values to a MAT file (optional)
save('groupedValues.mat', 'groupedValues');


%% To access the data 
keys = groupedValues.keys;
disp(keys)

value = groupedValues('Default'); % Replace 'key_name' with the actual key
disp(value);



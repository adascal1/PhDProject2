clc; clear all; close all; 

parcelPath = '/data/mica1/01_programs/micapipe-v0.2.0/parcellations';

%% Schaefer 100: 
schaefer100 = table2array(readtable(fullfile(parcelPath, 'schaefer-100_conte69.csv')));
    VISUAL      = [1:9, 51:58];
    SOMMOT      = [10:15, 59:66];
    DORSATT     = [16:23, 67:73];
    SALVENTATTN = [24:30, 74:78];
    LIMBIC      = [31:33, 79:80];
    CONT        = [34:37, 81:89];
    DEFAULT     = [38:50, 90:100];

    VISUAL100      = find(ismember(schaefer100, VISUAL));
    SOMMOT100      = find(ismember(schaefer100, SOMMOT));
    DORSATT100     = find(ismember(schaefer100, DORSATT));
    SALVENTATTN100 = find(ismember(schaefer100, SALVENTATTN));
    LIMBIC100      = find(ismember(schaefer100, LIMBIC));
    CONT100        = find(ismember(schaefer100, CONT));
    DEFAULT100     = find(ismember(schaefer100, DEFAULT));



% yeo100 = zeros(size(schaefer100));
% yeo100(VISUAL100)       = 1;
% yeo100(SOMMOT100)       = 2;
% yeo100(DORSATT100)      = 3 ;
% yeo100(SALVENTATTN100)  = 4;
% yeo100(LIMBIC100)       = 5;
% yeo100(CONT100)         = 6;
% yeo100(DEFAULT100)      = 7;




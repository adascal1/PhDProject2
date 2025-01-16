clc; clear all; close all; 

parcellation = strsplit(input("What threshold \n" , "s"));
    thresh   = str2double(parcellation{1});

 info = readtable(sprintf('/mfip/mfip1/arielle/PhDProject2/results/MICs/degree/degree_schaefer100_cohendResults_thresh-%f.csv', thresh));


degree100vs300 = find(info.significance == 1);
degree100vs300 = num2str(degree100vs300', '%g,')


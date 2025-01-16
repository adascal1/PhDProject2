clc; clear all; close all; 

info = readtable('/mfip/mfip1/arielle/PhDProject2/results/MICs/consistency/intersectionOverUnion_interParcellation.csv');

parcellation = strsplit(input("What parcellation schemes? (i.e.: 100 300) \n" , "s"));
    parcelNum1   = str2double(parcellation{1});
    parcelNum2   = str2double(parcellation{2});
    

degree100vs300 = find(info.ParcelNumber1 == parcelNum1 & info.ParcelNumber2 == parcelNum2 & strcmp(info.metric, 'eigenvectorCentrality'));
degree100vs300 = info(degree100vs300, :); 
degree100vs300 = find(degree100vs300.intersectionOverUnion > 0.5);
degree100vs300 = num2str(degree100vs300', '%g,')


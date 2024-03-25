% NMR Data Analysis Script
%
% This script reads NMR data from an Excel file \collated data\NMR_data.xlsx, 
% performs molecular mixing model calculations, and generates fractional
% compositions of different five molecular classes. It also saves the results to a 
% new Excel file (NMR_MMM.xlsx) and plots the fractional compositions.
% Usage:
%   1. Make sure the Excel file containing NMR data is located in the
%      directory specified by 'fname'.
%   2. Run the script in MATLAB.
%
% Author: Arjun Chakrawal
% Affiliation: Stockhom Univeristy
% Contact: arjun.chakrawal@natgeo.su.se
% For more information or questions, please contact the author.
%

clear all;
close all;clc

fname= "..\collated data\NMR_data.xlsx";

rawdata=readtable(fname,Sheet="all_data");
sz= size(rawdata);
rawdata.carbohydrate_MMM=rawdata.carbohydrate;
rawdata.protein_MMM=rawdata.protein;
rawdata.lignin_MMM=rawdata.lignin;
rawdata.lipid_MMM=rawdata.lipid;
rawdata.carbonyl_MMM=rawdata.carbonyl;
rawdata.char_MMM=rawdata.char;
rawdata.C=zeros(sz(1),1);rawdata.N=zeros(sz(1),1);rawdata.H=zeros(sz(1),1);rawdata.O=zeros(sz(1),1);
rawdata.molecularFormula = " " + zeros(sz(1),1);
rawdata.Cox = zeros(sz(1),1);
rawdata.r2 = zeros(sz(1),1);
rawdata.rmse = zeros(sz(1),1);
id = isnan(rawdata.C_N);

%% test 
temp = rawdata(strcmp(rawdata.Study,'Xu et al 2017 SBB'),:);

clc
i=1;
NMR_data=temp(i,:);
NMR_data.C_N
NMR_data.system
NMR_data.source
NMR_data.Csource
[frac, CNHO,  rmse, r_squared, new_mat]=molecular_mixing_model(NMR_data);
%% 

sz= size(rawdata);
for i =1:sz(1)
    i
    NMR_data=rawdata(i,:);
    if(isnan(NMR_data.carbohydrate) && ~isnan(NMR_data.AALKYL0_45Ppm))
        [frac,CNHO,  rmse, r_squared]=molecular_mixing_model(NMR_data);
        if(numel(frac)==4)
            rawdata.carbohydrate_MMM(i)=frac(1);
            rawdata.protein_MMM(i)=frac(2);
            rawdata.lipid_MMM(i)=frac(3);
            rawdata.carbonyl_MMM(i)=frac(4);
        else
            rawdata.carbohydrate_MMM(i)=frac(1);
            rawdata.protein_MMM(i)=frac(2);
            rawdata.lignin_MMM(i)=frac(3);
            rawdata.lipid_MMM(i)=frac(4);
            rawdata.carbonyl_MMM(i)=frac(5);
        end
        if(numel(frac)==6)
            rawdata.char_MMM(i)=frac(6);
        end
        rawdata.r2(i,:)= r_squared;
        rawdata.rmse(i,:) = rmse;

        rawdata.C(i,:)=CNHO.C;rawdata.N(i,:)=CNHO.N;
        rawdata.H(i,:)=CNHO.H;rawdata.O(i,:)=CNHO.O;
        rawdata.molecularFormula(i,:) =CNHO.molecularFormula;
        rawdata.Cox(i,:) = CNHO.Cox;
    end
end
temp = sum([rawdata.carbohydrate_MMM,rawdata.protein_MMM,rawdata.lignin_MMM,...
    rawdata.lipid_MMM,rawdata.carbonyl_MMM],2,"omitnan");
delete("../collated data/NMR_MMM.xlsx")
writetable(rawdata,"../collated data/NMR_MMM.xlsx")

fractions = [rawdata.carbohydrate_MMM,rawdata.protein_MMM,rawdata.lignin_MMM,...
    rawdata.lipid_MMM,rawdata.carbonyl_MMM, rawdata.char_MMM];
fractions=fractions(~isnan(fractions(:,1)),:);

fig=figure;fig.Color='w'; 
bar(fractions,'stacked', 'BarWidth',1)
ylim([0,1.1])
lh=legend("Carbohydrate","Protein","Lignin","Lipid","Carbonyl","Char");
lh.Location="eastoutside";lh.Box='off';lh.FontSize=12;





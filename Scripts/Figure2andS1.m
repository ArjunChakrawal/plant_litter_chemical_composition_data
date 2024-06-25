clear all; 
close all; clc
%%
rawdata=readtable("../collated data/NMR_MMM.xlsx");

sz=size(rawdata);
rawdata.ID = (1:sz(1))';
% Converting mg/g to grams per glitter
rawdata.Extract1_NPE = rawdata.Extract1_NPE *0.001; % g/glitter
rawdata.Extract2_WSE = rawdata.Extract2_WSE *0.001; % g/glitter
rawdata.Extract3_ACID = rawdata.Extract3_ACID *0.001; % g/glitter
rawdata.Extract4_AUR = rawdata.Extract4_AUR *0.001; % g/glitter

rawdata.AS = rawdata.Extract3_ACID; % g/glitter
rawdata.AIS = rawdata.Extract4_AUR; % g/glitter
disp(" Number of studies with Plant data " + length(unique(rawdata.Study)))

rawdata.labile =  sum(rawdata{:, {'Extract1_NPE', 'Extract2_WSE'}},2,"omitnan");
rawdata.labile(rawdata.labile==0)=nan; 

% filter rawdata for terrestrial system
terr_data = rawdata(strcmp(rawdata.system, 'terrestrial'),:);
% filter terrestrial system data only for plant material
plant_data0=terr_data(strcmp(terr_data.source,"plantC"),:);
plant_data=plant_data0;

disp(" Number of studies with Plant data " + length(unique(plant_data.Study)))

plant_data.litter_info = strcat(plant_data.type_Site_," ", plant_data.Species);
% filter plant_data where we have NMR observations
ind = any(~isnan(plant_data{:, 'carbohydrate_MMM'}), 2);
plant_data=plant_data(ind,:);
disp(" Number of studies with Plant data" + length(unique(plant_data.Study)))

% filter plant_data where we at least one proximate fraction observations
ind = any(~isnan(plant_data{:, {'labile','AS', 'AIS'}}), 2);
sum(ind)
plant_data=plant_data(ind,:);
disp(" Number of studies with Plant data" + length(unique(plant_data.Study)))

% remove rows where both AS and AIS are not available 
ind = all(~isnan(plant_data{:, {'AS', 'AIS'}}), 2);
sum(ind)
plant_data=plant_data(ind,:);
disp(" Number of studies with Plant data" + length(unique(plant_data.Study)))


% renomarmalize the sum PA fractions to 1 where sum >0.9
% create a new table with these corrected observations PD_PAsum0_9
S= sum(plant_data{:, {'labile','AS', 'AIS'}},2);
S_id= S>0.9;sum(S_id);
PD_PAsum0_9 = plant_data(S_id,:);

S= sum(PD_PAsum0_9{:, {'labile','AS', 'AIS'}},2);
PD_PAsum0_9.labile = PD_PAsum0_9.labile./S;
PD_PAsum0_9.AS = PD_PAsum0_9.AS./S;
PD_PAsum0_9.AIS = PD_PAsum0_9.AIS./S;
sum(PD_PAsum0_9{:, {'labile','AS', 'AIS'}},2)
disp(" Number of studies with Plant data" + length(unique(PD_PAsum0_9.Study)))

% remove PD_PAsum0_9 rows from orignal table 
plant_data(S_id,:)=[];

% for all rows where labile PA is not available, calculate it using
% 1-AS-AIS, and create a new table that contains all these rows PD_PA_newlabile
ind = all(~isnan(plant_data{:, {'AS', 'AIS'}}), 2);
sum(ind)

PD_PA_newlabile = plant_data(ind,:);
PD_PA_newlabile.labile= 1-PD_PA_newlabile.AS-...
    PD_PA_newlabile.AIS;
S=sum(PD_PA_newlabile{:, {'AS', 'AIS'}},2,"omitnan");
S_id= S<0.2;sum(S_id)
PD_PA_newlabile(S_id,:)=[];
disp(" Number of studies with Plant data" + length(unique(PD_PA_newlabile.Study)))

% remove PD_PA_newlabile rows from orignal table 
plant_data(ind,:)=[];

% get rid of row where we dont have AIS PA fraction
ind=isnan(plant_data.AIS); plant_data(ind,:)=[];
% now the remaining rows in the plant_data are either nan for AS or AIS


% combine 
PA_plant = [PD_PAsum0_9;PD_PA_newlabile];
PA_plant.sum=sum(PA_plant{:, {'labile','AS', 'AIS'}},2);
disp(" Number of studies with Plant data" + length(unique(PA_plant.Study)))

% remove PA or Biochar row from PA_plant
Csource_filter = ~ismember(PA_plant.Csource, {'PA', 'biochar'});
PA_plant = PA_plant(Csource_filter, :);
disp(" Number of studies with Plant data" + length(unique(PA_plant.Study)))

PA_plant= PA_plant(~isnan(PA_plant.carbohydrate_MMM),:);
PA_plant= PA_plant(~isnan(PA_plant.lignin_MMM),:);
% remove row with unrealistic %C from PA_plant
id=PA_plant.C_conc>=0.3;
PA_plant=PA_plant(id,:);
disp(" Number of studies with Plant data" + length(unique(PA_plant.Study)))

%%
% Carbohydrate	Protein	Lignin	Lipid	Carbonyl	Char
CNHO_terrestrial = [1, 1, 1, 1, 1, 1; ...
    0, 0.27, 0, 0, 0, 0; ...
    1.67, 1.10, 1.24, 1.94, 1, 0.45; ...
    0.83, 0.16, 0.43, 0.24, 2, 0.41];
gCNMR= 12./(CNHO_terrestrial'*[12,14,1,16]');

PA_plant.Carbohydrates=PA_plant.carbohydrate_MMM.*PA_plant.C_conc.*0.001./gCNMR(1); % gC/gC litter * gC litter/glitter /gC/gMolecule    
PA_plant.Proteins=PA_plant.protein_MMM.*PA_plant.C_conc.*0.001./gCNMR(2); % g/glitter
PA_plant.Lignins=PA_plant.lignin_MMM.*PA_plant.C_conc.*0.001./gCNMR(3); % g/glitter
PA_plant.Lipids=PA_plant.lipid_MMM.*PA_plant.C_conc.*0.001./gCNMR(4); % g/glitter
PA_plant.Carbonyls=PA_plant.carbonyl_MMM.*PA_plant.C_conc.*0.001./gCNMR(5); % g/glitter

S= sum(PA_plant{:, {'Carbohydrates', 'Proteins','Lignins', 'Lipids', 'Carbonyls'}},2);
S= sum(PA_plant{:, {'carbohydrate_MMM', 'protein_MMM','lignin_MMM', 'lipid_MMM', 'carbonyl_MMM'}},2);

col2={'ID','Latitude','Longitude','Csource','Study','litter_info',...
    'timeDay','C_conc','C_N','Extract1_NPE', 'Extract2_WSE','Extract3_ACID','Extract4_AUR',...
    'labile','AS', 'AIS','AALKYL0_45Ppm','BMETHOX45_60Ppm','CO_ALKYL60_95Ppm','DDI_O_ALK95_110ppm',...
    'EAROM110_145Ppm','FPHEN145_165Ppm','GCARBOX165_210Ppm','Carbohydrates', 'Proteins',...
    'Lignins', 'Lipids', 'Carbonyls','Cox','r2','rmse'};

cordata = PA_plant(:,col2);
id = find(strcmp(cordata.Csource,'cropresidue') | strcmp(cordata.Csource,'grass'));
for i = 1:length(id)
    cordata.Csource{id(i)} = 'grass+herbs';
end
unique(cordata.Csource)

sz=size(cordata);
cordata.ID = (1:sz(1))';
study  = unique(cordata.Study);
delete("../collated data/corrdata.xlsx");
writetable(cordata,'../collated data/corrdata.xlsx');

%%


figure;
boxplot([cordata.r2, cordata.rmse], 'Labels', {'r2', 'rmse'});
title("molecular mixing model performance")

vars={'labile','AS','AIS', ...
    'Carbohydrates', 'Proteins','Lignins', 'Lipids', 'Carbonyls'};

[fig, ~,~]=paired_correlation(cordata,vars,vars, 0.05,'', ...
    20,'MarkerFaceColor',[255, 145, 0]/255,'MarkerFaceAlpha',0.1, ...
    'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1);
fig.Position = [200 60  888   794];
exportgraphics(fig,'../results/FigureS2.png',Resolution=300)


x_vars =  vars(1:3);
y_vars = vars(4:end);
LC=summer(length(cordata.AS));
[fig, ax,t]=paired_correlation(cordata,x_vars,y_vars, 0.05,'', ...
    20,'MarkerFaceColor',LC(1,:),'MarkerFaceAlpha',0.2, ...
    'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1);
fig.Position(3:4) = [700  700.60];
xlab = {'Labile (nonpolar and water soluble)','Acid soluble (AS)','Acid insoluble (AIS)'};

for i =1:3
    xlabel(ax(5,i),xlab{i},'fontsize',11)
end
xlabel(t,'Fractions of extractives from proximate analysis (g/g litter)')
ylabel(t,'Fractions of organic compounds estimated from NMR (g/g litter)')

exportgraphics(fig,'../results/Figure2.png',Resolution=300)


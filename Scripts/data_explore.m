

clear all;
close all;clc
rawdata=readtable("../collated data/NMR_MMM.xlsx", Sheet="Sheet1");
unique(rawdata.Study)
%% geolocation map
systValues = unique(rawdata.system(~isnan(rawdata.Longitude)));
srcdata=rawdata(strcmp(rawdata.source, "plantC"),:);
lat=srcdata.Latitude(~isnan(srcdata.Latitude));
long=srcdata.Longitude(~isnan(srcdata.Longitude));
figure;geoscatter(lat, long, 20,"filled"); hold on
geobasemap bluegreen;set(gcf,'color','w')
exportgraphics(gcf, '../results/plantC_latlong.png', Resolution=300)


for i =1:numel(systValues)
    fig=figure; tiledlayout(4,2,"TileSpacing","compact","Padding","compact");
    fig.Position=[200 200 1155  764];fig.Color='w';
    data=rawdata(strcmp(rawdata.system, systValues{i}),:);
    SrcValues = unique(data.source);
    for j =1:numel(SrcValues)
        srcdata=data(strcmp(data.source, SrcValues{j}),:);
        lat=srcdata.Latitude(~isnan(srcdata.Latitude));
        long=srcdata.Longitude(~isnan(srcdata.Longitude));
        nexttile;geoscatter(lat, long, 20,"filled"); hold on
        title(SrcValues{j}+", n="+size(srcdata,1))
        geolimits([-50 80],[-180 180])
    end
end
exportgraphics(fig, '../results/datamap.png', Resolution=300)
%% barplot all fractions
fractions = [rawdata.carbohydrate_MMM,rawdata.protein_MMM,rawdata.lignin_MMM, rawdata.lipid_MMM,rawdata.carbonyl_MMM, rawdata.char_MMM];
fractions=fractions(~isnan(fractions(:,1)),:);

fig=figure;fig.Color='w';fig.Position=[183 557 1439 408];
bar(fractions,'stacked','BarWidth', 1)
ylim([0,1])
lh=legend("Carbohydrate","Protein","Lignin","Lipid","Carbonyl","Char");
lh.Location="eastoutside";lh.Box='off';lh.FontSize=12;
exportgraphics(fig, '../results/barplot_fractions.png', Resolution=300)
%% plot molecular fractions
terr_data = rawdata(strcmp(rawdata.system, 'terrestrial'),:);
grpstats(terr_data,"source")
grpstats(terr_data,"Csource")
srcVal=unique(terr_data.source);
srcVal(strcmp(srcVal, 'MicrobialC'))=[];

fig=figure;fig.Color='w';fig.Position=[200         122        1080         637];
tiledlayout('flow',TileSpacing='tight',Padding='tight')
for i =1:numel(srcVal)
    data=terr_data(strcmp(terr_data.source, srcVal{i}),:);
    Csrc=unique(data.Csource);
    for j =1:numel(Csrc)
        C_data =data(strcmp(data.Csource, Csrc{j}),:);
        fractions = [C_data.carbohydrate_MMM,C_data.protein_MMM,C_data.lignin_MMM,C_data.lipid_MMM,C_data.carbonyl_MMM,C_data.char_MMM];
        fractions=fractions(~isnan(fractions(:,1)),:);
        fractions=fractions(:,~isnan(fractions(1,:)));
        if(size(fractions,1)==1)
            nexttile;bar(1,fractions,'stacked','BarWidth', 1,EdgeColor='none');
        else
            nexttile;bar(fractions,'stacked','BarWidth', 1,EdgeColor='none');
        end
        title(Csrc{j}+" (n="+size(C_data,1)+")"); axis tight
    end
end
lh=legend("Carbohydrate","Protein","Lignin","Lipid","Carbonyl","Char");
lh.Location="eastoutside";lh.Box='off';lh.FontSize=11;
lh.Position=[ 0.4334    0.0255    0.1370    0.2080];
allaxes = findall(fig, 'type', 'axes');
set(allaxes,box='off', yticklabels=[]);
exportgraphics(fig, '../results/terrestrial_fractions.png', Resolution=300)



aq_data = rawdata(strcmp(rawdata.system, 'aquatic'),:);
grpstats(aq_data,"source")
grpstats(aq_data,"Csource")
srcVal=unique(aq_data.source);
srcVal(strcmp(srcVal, 'MicrobialC'))=[];

fig=figure;fig.Color='w';fig.Position=[200         122        876   424];
tiledlayout('flow',TileSpacing='tight',Padding='tight')
for i =1:numel(srcVal)
    data=aq_data(strcmp(aq_data.source, srcVal{i}),:);
    Csrc=unique(data.Csource);
    for j =1:numel(Csrc)
        C_data =data(strcmp(data.Csource, Csrc{j}),:);
        fractions = [C_data.carbohydrate_MMM,C_data.protein_MMM,C_data.lignin_MMM,C_data.lipid_MMM,C_data.carbonyl_MMM,C_data.char_MMM];
        fractions=fractions(~isnan(fractions(:,1)),:);
        fractions=fractions(:,~isnan(fractions(1,:)));
        nexttile;bar(fractions,'stacked','BarWidth', 1,EdgeColor='none');
        title(Csrc{j}+" (n="+size(C_data,1)+")"); axis tight
    end
end
lh=legend("Carbohydrate","Protein","Lignin","Lipid","Carbonyl","Char");
lh.Location="eastoutside";lh.Box='off';lh.FontSize=12;
exportgraphics(fig, '../results/aquatic_fractions.png', Resolution=300)


fig=figure;fig.Color='w';fig.Position=[200         122        876   424];
tiledlayout('flow',TileSpacing='tight',Padding='tight')
microbe_data = rawdata(strcmp(rawdata.source, 'MicrobialC'),:);
Csrc=unique(microbe_data.Csource);
for j =1:numel(Csrc)
    C_data =microbe_data(strcmp(microbe_data.Csource, Csrc{j}),:);
    fractions = [C_data.carbohydrate_MMM,C_data.protein_MMM,C_data.lignin_MMM,C_data.lipid_MMM,C_data.carbonyl_MMM,C_data.char_MMM];
    fractions=fractions(~isnan(fractions(:,1)),:);
    fractions=fractions(:,~isnan(fractions(1,:)));
    if(size(fractions,1)==1)
        nexttile;bar(1,fractions,'stacked','BarWidth', 1,EdgeColor='none');
    else
        nexttile;bar(fractions,'stacked','BarWidth', 1,EdgeColor='none');
    end
    title(Csrc{j}+" (n="+size(C_data,1)+")"); axis tight
end
lh=legend("Carbohydrate","Protein","Lipid","Carbonyl");
lh.Location="eastoutside";lh.Box='off';lh.FontSize=12;
exportgraphics(fig, '../results/microbial_fractions.png', Resolution=300)

%% plot oxidation state
close all
fig=figure;fig.Color='w';
plotSwarmChart(rawdata,"Cox","system",'o','filled','MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.3)
ylabel("Oxidation state"); title('all data')
exportgraphics(fig, '../results/all_DR1.png', Resolution=300)


fig=figure;fig.Color='w';
plotSwarmChart(rawdata,"Cox","source",'o','filled','MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.3)
ylabel("Oxidation state"); title('all data')
exportgraphics(fig, '../results/all_DR2.png', Resolution=300)

rawdata.Csource=categorical(rawdata.Csource);
fig=figure;fig.Color='w';fig.Position=[200 200 1000 344];
fresh= rawdata(rawdata.timeDay==0,:);
decom= rawdata(rawdata.timeDay~=0,:);
swarmchart(fresh.Csource,fresh.Cox,'o','filled','MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.3); hold on
swarmchart(decom.Csource,decom.Cox,'o','filled','MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.15);
lh=legend("Fresh","Decomposed"); lh.Box='off';lh.FontSize=12;
set(gca, 'TickLabelInterpreter', 'none', 'FontSize',12);
ylabel("Oxidation state"); title('all data'); grid on 
exportgraphics(fig, '../results/all_DR3.png', Resolution=300)


fig=figure;fig.Color='w';fig.Position=[200   200   680   455];
terr_data = rawdata(strcmp(rawdata.system, 'terrestrial'),:);
plotSwarmChart(terr_data,"Cox","source",'o','filled','MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.3)
set(gca, 'TickLabelInterpreter', 'none', 'FontSize',12);
ylabel("Oxidation state"); title('terrestrial')
exportgraphics(fig, '../results/terr_degree_reduc.png', Resolution=300)

% -------------------------------------------------------------------------

terr_data = rawdata(strcmp(rawdata.system, 'terrestrial'),:);
plant_data = terr_data(strcmp(terr_data.source, 'plantC'),:);
grpstats(plant_data,"Study")

plant_data(strcmp(plant_data.Csource, 'PA'),:)=[];
plant_data.Csource=categorical(plant_data.Csource);


fig=figure;fig.Color='w';fig.Position=[200   200   680   455];
tiledlayout('flow',TileSpacing='tight',Padding='tight')
fresh= plant_data(plant_data.timeDay==0,:);
decom= plant_data(plant_data.timeDay~=0,:);
swarmchart(fresh.Csource,fresh.Cox,'o','filled','MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.15); hold on
swarmchart(decom.Csource,decom.Cox,'o','filled','MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.3);
ylabel("Oxidation state");title('terrestrial plantC')
lh=legend("Fresh","Decomposed"); lh.Box='off';lh.FontSize=12;
set(gca, 'TickLabelInterpreter', 'none', 'FontSize',12);
exportgraphics(fig, '../results/terr_plant_degree_reduc.png', Resolution=300)

% -------------------------------------------------------------------------
aq_data = rawdata(strcmp(rawdata.system, 'aquatic'),:);
fig=figure;fig.Color='w';fig.Position=[200   200   680   455];
plotSwarmChart(aq_data,"Cox","source",'o','filled','MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.3)
set(gca, 'TickLabelInterpreter', 'none', 'FontSize',12);
ylabel("Oxidation state");title('aquatic')
exportgraphics(fig, '../results/aq_degree_reduc.png', Resolution=300)


aq_data = rawdata(strcmp(rawdata.system, 'aquatic'),:);
plant_data = aq_data(strcmp(aq_data.source, 'plantC'),:);
plant_data.Csource=categorical(plant_data.Csource);
fig=figure;fig.Color='w';fig.Position=[200   200   680   455];
tiledlayout('flow',TileSpacing='tight',Padding='tight')
fresh= plant_data(plant_data.timeDay==0,:);
decom= plant_data(plant_data.timeDay~=0,:);
swarmchart(fresh.Csource,fresh.Cox,'o','filled','MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.15); hold on
swarmchart(decom.Csource,decom.Cox,'o','filled','MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.3);
ylabel("Oxidation state");title('aquatic')
lh=legend("Fresh","Decomposed"); lh.Box='off';lh.FontSize=12;
set(gca, 'TickLabelInterpreter', 'none', 'FontSize',12);
exportgraphics(fig, '../results/aq_plant_degree_reduc.png', Resolution=300)


%%











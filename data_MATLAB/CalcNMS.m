%Plant project with Jennifer Bowen's group at Northeastern
%Data are Lumos, negative ion mode, untargeted analysis
%Calculate the NMS ordinations to provide a quick look at the pattern
%Krista Longnecker 2/21/2011; updated 10/11/2023 to use for plant project
%
%notes on functions...this code will require a few things that I am not 
%putting onto GitHub:
% 1. The MATLAB statistics toolox
% 2. The Fathom toolbox which is available at Mathworks file exchange:
% https://www.mathworks.com/matlabcentral/fileexchange/68518-fathom-toolbox-for-matlab
% 3. I use cbrewer for colors (also available at the Mathworks file
% exchange); though can simply change to a different option for colors
%https://www.mathworks.com/matlabcentral/fileexchange/58350-cbrewer2?s_tid=srchtitle
% 4. title_up.m which I will upload to GitHub, but is a simple way to put a
% title centered on the plots. Perhaps MATLAB has corrected this, but it
% used to be harder to do this nicely.

clear all
close all
load NEplants_neg_aligned.2023.12.06.mat  

%set up some trimming to allow only a subset of the data to be analyzed
if 0
    %use everything (probably not desired)
    Intensity = EICdata;
    ks = [1:size(EICdata,2)];
elseif 0
    %plot unknowns & pooled
    s = strcmp(sInfo.sample,'Unknown');
    sp = strcmp(sInfo.sample,'pooled');
    ks = find(s==1 | sp==1);
    Intensity = EICdata(:,ks);
    sInfo = sInfo(ks,:);
    clear s sp ks
elseif 1
    %only plot the Unknowns
    s = strcmp(sInfo.sample,'Unknown');
    ks = find(s==1);
    Intensity = EICdata(:,ks);
    sInfo = sInfo(ks,:);
    clear s sp ks
    
end
clear fileName EICdata

%change any NaNs to zero
i = isnan(Intensity);
ki = find(i==1);
Intensity(ki) = 0;

UseIntensity=2; %decide which option to use for peak heights
switch UseIntensity
    case 1  % - don't convert except to presence/absence
        fm=f_normal(Intensity,'01');
    case 2
        %how to transform, perhaps just sqrt
        fm=f_transform(Intensity,1);
end
      
Transform=1;    %decide on data transformation and distance measure
switch Transform
    case 1
        %just do the Bray-Curtis distance measure
        xdist=f_braycurtis(fm);
        useT=(['P/A Bray Curtis']);
    case 2
        %first relativize the data by SU totals, then do Bray-Curtis
        t=f_transform(fm,6);
        xdist=f_braycurtis(t);
        useT=([' relative, P/A Bray Curtis']);
end
clear t fm*

if 1 %set to 0 to skip the NMS
    opts = statset('display', 'final', 'maxiter', 200,'tolFun',1e-6);
    [Y , stress] =mdscale(xdist,2,'start','random','options',opts,'criterion','stress'); %assuming two dimensions
    clear opts

    %calculate the stats for the NMS: total variabilty and for each axis
    r2=(f_mantel(xdist,f_euclid(Y'),1).r).^2
    r2axis1=(f_mantel(xdist,f_euclid(Y(:,1)'),1).r).^2
    r2axis2=(f_mantel(xdist,f_euclid(Y(:,2)'),1).r).^2

    figure
    subplot(2,4,[1 2 5 6])
    plot(Y(:,1),Y(:,2),'k.')
    hold on
    %label each point on the NMS with the GenotypeID from Matt
    for i=1:length(xdist)
        h=text(Y(i,1)+0.01 , Y(i,2) , sInfo.GenotypeID(i)) ;
        set(h,'FontSize',10,'Color','k')
        clear h
    end;
    clear i
    title('Genotype ID')

    %shuffle colors manually...
    if 0
        cmap = cbrewer('qual','Dark2',8);
        cmap = cmap([1 5 2 3 4 6:8],:);
    elseif 1
        cmap = cbrewer('qual','Set1',9);
        cmap = cmap([1 8 2 3 4:7 9],:);
    end
    
    %set up variables for the markers
    ms=10;
    um = {'v','<','p','h'};
    
    subplot(243)
    h = gscatter(Y(:,1),Y(:,2),{sInfo.Site sInfo.Ecotype},parula(5),[],30);
    for ah = 1:length(h)
        set(h(ah),'markerfacecolor',cmap(ah,:),'markeredgecolor',cmap(ah,:),...
            'marker',um{ah},'markersize',ms)
    end
    clear ah
    title('site x ecotype')
    
    subplot(247)
    h = gscatter(Y(:,1),Y(:,2),sInfo.Site,parula(5),[],30);
    cmap = cbrewer('qual','Set1',9);
    cmap = cmap([1 2 3 4:9],:);    
    um = {'v','p','h'};    
    for ah = 1:length(h)
        set(h(ah),'markerfacecolor',cmap(ah,:),'markeredgecolor',cmap(ah,:),...
            'marker',um{ah},'markersize',ms)
    end
    clear ah    
    title('sites')
    
    subplot(248)
    h = gscatter(Y(:,1),Y(:,2),sInfo.Ecotype,parula(5),[],30);
    cmap = cbrewer('qual','Dark2',8);
    um = {'s','*','h'};    
    for ah = 1:length(h)
        set(h(ah),'markerfacecolor',cmap(ah,:),'markeredgecolor',cmap(ah,:),...
            'marker',um{ah},'markersize',ms)
    end
    clear ah    
    title('ecotype')
    
    subplot(244)
    h = gscatter(Y(:,1),Y(:,2),sInfo.dead,parula(5),[],30);
    cmap = cbrewer('div','BrBG',3);
    cmap = cmap([1 3 2],:);    
    um = {'+','s'};   
    for ah = 1:length(h)
        set(h(ah),'markeredgecolor',cmap(ah,:),...
            'marker',um{ah},'markersize',ms,'lineWidth',3)
    end
    clear ah    
    title('Colors are notes about potentially dead plants')   

    set(gcf,'paperpositionmode','auto','position',[-1764 72 1644 927])

end

title_up('Bowen samples, round 2, October 2023')
set(gcf,'paperpositionmode','auto')
saveas(gcf,'Bowen_plants2_NMSplots.2023.12.06.pdf','pdf')


if 0
    %use this to plot a cluster diagram
    figure
    Z=linkage(squareform(xdist),'average');

    StringLabels = sInfo.Bowen_name;
    dendrogram(Z,'Labels',StringLabels,'orientation','left');
    title_up('Bowen samples, round 2 October 2023')

    set(gcf,'position',[-159 69 1188 652])
    set(gcf,'paperpositionmode','auto')
    saveas(gcf,'Bowen_plants2_clusterDiagram.pdf','pdf')
end

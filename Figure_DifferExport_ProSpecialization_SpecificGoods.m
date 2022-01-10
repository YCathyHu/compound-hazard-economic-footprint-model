%% Figures for differexport on specific goods
% 30%-24 COVID control over all regions
% Multi-scale flooding in region C for 2 weeks
% Trade scenarios: 50% export restriction on MANR-C；MANR-C, MANR-B
% Specialized production: MANR-C, MANR-B
% 3*3 subplots reporting GDP variations in four regions with different
% combinations of flood and export restriction scenarios
% close all
clear
clc
%% Figure drawing
load FigureData_differexportsp.mat
load FigureData_basefrtrade_sp.mat

dx = 1:T;

linespec = ["-","-","-","-"];
colorspec = [...
    182, 115, 101;...
    100, 181, 205;...
    204, 185, 116;...
    184, 135, 195;...
    164, 160, 155]./255;
colorbg = [...
    242, 242, 242;...
    217, 217, 217;...
    191, 191, 191;...
    166, 166, 166;...
    128, 128, 128]./255;
linewidthspec = [3,1,0.5];

t = ["Small flood in C", "Medium flood in C","Large flood in C"];

figure(1);
clf;

set(gcf,'Units','centimeter','Position',[1 0.8 15.92 9.84]); % set figure size

f = tiledlayout(3,3);
f.TileSpacing = 'compact';
f.Padding = 'compact';

% Free trade scenarios
for i = 1:size(Flood_set_int,1)
    nexttile
    dy = 100.*[zeros(1,T);VA_AT3(i,:);VA_OT3(i,:);VA_CT3(i,:);VA_DT3(i,:)];
    p1 = plot(dx,dy,'LineWidth', linewidthspec(2));
    for j = 1:R+1
        if j ==1
            p1(j).Color = colorspec(end,:);
            p1(j).LineStyle = '--';
        else
            p1(j).Color = colorspec(j-1,:);
            p1(j).LineStyle = linespec(1);
        end
    end
    hold on
    ylim([-80,10])
    xlim([0,T])
    set(gca,'YTick',-60:30:0,'XTick',0:20:80);
    ylb = get(gca,'YTickLabel');
    n = length(ylb);
    pp = '%';
    new_ylb = cell(n,1);
    for k = 1:n
        new_ylb{k}= [ylb{k},pp];
    end
    set(gca,'yticklabel',new_ylb,'Layer','top','FontName','Helvetica','Fontsize',8,'LineWidth',linewidthspec(3),'Color',[colorbg(2,:),0.3]);
    box off
    title(t(i))
end

% Export restriction scenarios
for sc = 1:size(Scenarios,1)
    for i = 1:size(Flood_set_int,1)
        nexttile
        dy = 100.*[zeros(1,T);VA_AT(i+(sc-1)*size(Flood_set_int,1),:);VA_OT(i+(sc-1)*size(Flood_set_int,1),:);VA_CT(i+(sc-1)*size(Flood_set_int,1),:);VA_DT(i+(sc-1)*size(Flood_set_int,1),:)];
        p1 = plot(dx,dy,'LineWidth', linewidthspec(2));
        for j = 1:R+1
            if j ==1
                p1(j).Color = colorspec(end,:);
                p1(j).LineStyle = '--';
            else
                p1(j).Color = colorspec(j-1,:);
                p1(j).LineStyle = linespec(1);
            end
        end
        hold on
        ylim([-80,10])
        xlim([0,T])
        set(gca,'YTick',-60:30:0,'XTick',0:20:80);
        ylb = get(gca,'YTickLabel');
        n = length(ylb);
        pp = '%';
        new_ylb = cell(n,1);
        for k = 1:n
            new_ylb{k}= [ylb{k},pp];
        end
        set(gca,'yticklabel',new_ylb,'Layer','top','FontName','Helvetica','Fontsize',8,'LineWidth',linewidthspec(3),'Color',[colorbg(2,:),0.3]);
        if sc == size(Scenarios,1)
            xlabel('Weeks');
        end
        if sc == 1 && i == 1
            ylabel('% change of regional GDP');
        end
        box off
    end
end
lgd = legend('Pre-disaster','Region A','Region B','Region C','Region D');
set(lgd,'NumColumns',1,'FontSize',7,'FontName','Helvetica','Box','off','Location','SouthEast');
%% End.

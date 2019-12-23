function plot_4_potential_dam_characteristics_portrait(dam_ids,breite,hoehe,schriftgroesse,visible)
% ######  fact sheet, page 2  ######
% plot the dam characteristics, basin characteristics and evaluation
% criteria for all potential dam heights for one or more dam sites
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019

global basins_selected

clear Legende p

if length(dam_ids)>1
    farbe = jet(length(basins_selected))*.8;
    for i_dam = dam_ids
        Legende{i_dam} = sprintf('Basin %d',i_dam);
    end
    Legende{end+1} = 'possible height';
    Legende{end+1} = 'selected height';
else
    i_dam = dam_ids;
    farbe(i_dam,:) = [0,0,0];
    Legende{1} = 'possible height';
    Legende{2} = 'selected height';
end


if visible == 0
    fig = figure('color',[1,1,1],'rend','painters','paperunits','centimeters','papersize',[breite,hoehe],...
        'paperpositionmode','manual','paperPosition',[0 0 breite hoehe],'position',[1 1 breite hoehe]*50,'visible','off');
else
    fig = figure('color',[1,1,1],'rend','painters','paperunits','centimeters','papersize',[breite,hoehe],...
            'paperpositionmode','manual','paperPosition',[0 0 breite hoehe],'position',[1 1 breite hoehe]*50);    
end


%% dam characteristics
set(0,'CurrentFigure',fig)

% dam volume
sp1 = subplot(3,4,1);
hold on
i_selected = find(basins_selected(dam_ids(1)).curve_potential_dam_basin(:,1)==basins_selected(dam_ids(1)).dam_height);
    i_possible = isnan(basins_selected(dam_ids(1)).curve_potential_dam_basin(:,8))==0;
p2 = plot(basins_selected(dam_ids(1)).curve_potential_dam_basin(i_possible,2),basins_selected(dam_ids(1)).curve_potential_dam_basin(i_possible,1),...
    '.','color',[.5,.5,.5]);
p3 = plot(basins_selected(dam_ids(1)).curve_potential_dam_basin(i_selected,2),basins_selected(dam_ids(1)).curve_potential_dam_basin(i_selected,1),...
    'o','color',[.5,.5,.5]);
for i_dam = dam_ids
    i_selected = find(basins_selected(i_dam).curve_potential_dam_basin(:,1)==basins_selected(i_dam).dam_height);
    i_possible = isnan(basins_selected(i_dam).curve_potential_dam_basin(:,8))==0;
    p(i_dam) = plot(basins_selected(i_dam).curve_potential_dam_basin(:,2),basins_selected(i_dam).curve_potential_dam_basin(:,1),...
        '-','color',farbe(i_dam,:));
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_possible,2),basins_selected(i_dam).curve_potential_dam_basin(i_possible,1),...
        '.','color',farbe(i_dam,:));
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_selected,2),basins_selected(i_dam).curve_potential_dam_basin(i_selected,1),...
        'o','color',farbe(i_dam,:));
end
xlabel('dam volume (m³)        ')
ylabel('dam height (m)')
grid on
grid minor


% dam axis length
sp2 = subplot(3,4,2);
hold on
for i_dam = dam_ids
    i_selected = find(basins_selected(i_dam).curve_potential_dam_basin(:,1)==basins_selected(i_dam).dam_height);
    i_possible = isnan(basins_selected(i_dam).curve_potential_dam_basin(:,8))==0;
    plot(basins_selected(i_dam).curve_potential_dam_basin(:,3),basins_selected(i_dam).curve_potential_dam_basin(:,1),...
        '-','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_possible,3),basins_selected(i_dam).curve_potential_dam_basin(i_possible,1),...
        '.','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_selected,3),basins_selected(i_dam).curve_potential_dam_basin(i_selected,1),...
        'o','color',farbe(i_dam,:))
end
xlabel('dam axis length (m)')
ylabel('dam height (m)')
grid on
grid minor


% real dam height
sp3 = subplot(3,4,3);
hold on
for i_dam = dam_ids
    i_selected = find(basins_selected(i_dam).curve_potential_dam_basin(:,1)==basins_selected(i_dam).dam_height);
    i_possible = isnan(basins_selected(i_dam).curve_potential_dam_basin(:,8))==0;
    plot(basins_selected(i_dam).curve_potential_dam_basin(:,9),basins_selected(i_dam).curve_potential_dam_basin(:,1),...
        '-','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_possible,9),basins_selected(i_dam).curve_potential_dam_basin(i_possible,1),...
        '.','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_selected,9),basins_selected(i_dam).curve_potential_dam_basin(i_selected,1),...
        'o','color',farbe(i_dam,:))
end
xlabel('dam height (real) (m)')
ylabel('dam height (m)')
grid on
grid minor


%% basin characteristics

% basin volume
sp5 = subplot(3,4,5);
hold on
for i_dam = dam_ids
    i_selected = find(basins_selected(i_dam).curve_potential_dam_basin(:,1)==basins_selected(i_dam).dam_height);
    i_possible = isnan(basins_selected(i_dam).curve_potential_dam_basin(:,8))==0;
    plot(basins_selected(i_dam).curve_potential_dam_basin(:,4),basins_selected(i_dam).curve_potential_dam_basin(:,1),...
        '-','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_possible,4),basins_selected(i_dam).curve_potential_dam_basin(i_possible,1),...
        '.','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_selected,4),basins_selected(i_dam).curve_potential_dam_basin(i_selected,1),...
        'o','color',farbe(i_dam,:))
end
xlabel('basin volume (m³)        ')
ylabel('dam height (m)')
grid on
grid minor

% basin area
sp6 = subplot(3,4,6);
hold on
for i_dam = dam_ids
    i_selected = find(basins_selected(i_dam).curve_potential_dam_basin(:,1)==basins_selected(i_dam).dam_height);
    i_possible = isnan(basins_selected(i_dam).curve_potential_dam_basin(:,8))==0;
    plot(basins_selected(i_dam).curve_potential_dam_basin(:,5),basins_selected(i_dam).curve_potential_dam_basin(:,1),...
        '-','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_possible,5),basins_selected(i_dam).curve_potential_dam_basin(i_possible,1),...
        '.','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_selected,5),basins_selected(i_dam).curve_potential_dam_basin(i_selected,1),...
        'o','color',farbe(i_dam,:))
end
xlabel('basin area (m²)')
ylabel('dam height (m)')
grid on
grid minor

% specific volume
sp7 = subplot(3,4,7);
hold on
for i_dam = dam_ids
    i_selected = find(basins_selected(i_dam).curve_potential_dam_basin(:,1)==basins_selected(i_dam).dam_height);
    i_possible = isnan(basins_selected(i_dam).curve_potential_dam_basin(:,8))==0;
    plot(basins_selected(i_dam).curve_potential_dam_basin(:,6),basins_selected(i_dam).curve_potential_dam_basin(:,1),...
        '-','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_possible,6),basins_selected(i_dam).curve_potential_dam_basin(i_possible,1),...
        '.','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_selected,6),basins_selected(i_dam).curve_potential_dam_basin(i_selected,1),...
        'o','color',farbe(i_dam,:))
end
xlabel('specific volume (mm)')
ylabel('dam height (m)')
grid on
grid minor


%% evaluation criteria

% objective function
sp8 = subplot(3,4,8);
hold on
for i_dam = dam_ids
    i_selected = find(basins_selected(i_dam).curve_potential_dam_basin(:,1)==basins_selected(i_dam).dam_height);
    i_possible = isnan(basins_selected(i_dam).curve_potential_dam_basin(:,8))==0;
    plot(basins_selected(i_dam).curve_potential_dam_basin(:,7),...
        basins_selected(i_dam).curve_potential_dam_basin(:,1),'-','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_possible,7),...
        basins_selected(i_dam).curve_potential_dam_basin(i_possible,1),'.','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_selected,7),...
        basins_selected(i_dam).curve_potential_dam_basin(i_selected,1),'o','color',farbe(i_dam,:))
end
xlabel('objective function (-)')
ylabel('dam height (m)')
grid on
grid minor

% dam volume per storage volume (criterion 1)
sp9 = subplot(3,4,9);
hold on
for i_dam = dam_ids
    i_selected = find(basins_selected(i_dam).curve_potential_dam_basin(:,1)==basins_selected(i_dam).dam_height);
    i_possible = isnan(basins_selected(i_dam).curve_potential_dam_basin(:,8))==0;
    plot(basins_selected(i_dam).curve_potential_dam_basin(:,10),...
        basins_selected(i_dam).curve_potential_dam_basin(:,1),'-','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_possible,10),...
        basins_selected(i_dam).curve_potential_dam_basin(i_possible,1),'.','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_selected,10),...
        basins_selected(i_dam).curve_potential_dam_basin(i_selected,1),'o','color',farbe(i_dam,:))
end
xlabel('dam volume/storage volume (-)')
ylabel('dam height (m)')
grid on
grid minor

% basin area per storage volume (criterion 2)
sp10 = subplot(3,4,10);
hold on
for i_dam = dam_ids
    i_selected = find(basins_selected(i_dam).curve_potential_dam_basin(:,1)==basins_selected(i_dam).dam_height);
    i_possible = isnan(basins_selected(i_dam).curve_potential_dam_basin(:,8))==0;
    plot(basins_selected(i_dam).curve_potential_dam_basin(:,11),...
        basins_selected(i_dam).curve_potential_dam_basin(:,1),'-','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_possible,11),...
        basins_selected(i_dam).curve_potential_dam_basin(i_possible,1),'.','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_selected,11),...
        basins_selected(i_dam).curve_potential_dam_basin(i_selected,1),'o','color',farbe(i_dam,:))
end
xlabel('basin area/storage volume (1/m)')
ylabel('dam height (m)')
grid on
grid minor

% share of well-suited cells (criterion 3)
sp11 = subplot(3,4,11);
hold on
for i_dam = dam_ids
    i_selected = find(basins_selected(i_dam).curve_potential_dam_basin(:,1)==basins_selected(i_dam).dam_height);
    i_possible = isnan(basins_selected(i_dam).curve_potential_dam_basin(:,8))==0;
    plot(basins_selected(i_dam).curve_potential_dam_basin(:,12),...
        basins_selected(i_dam).curve_potential_dam_basin(:,1),'-','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_possible,12),...
        basins_selected(i_dam).curve_potential_dam_basin(i_possible,1),'.','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_selected,12),...
        basins_selected(i_dam).curve_potential_dam_basin(i_selected,1),'o','color',farbe(i_dam,:))
end
xlabel('share of well-suited cells (-)')
ylabel('dam height (m)')
grid on
grid minor

% share of not-suited cells (criterion 4)
sp12 = subplot(3,4,12);
hold on
for i_dam = dam_ids
    i_selected = find(basins_selected(i_dam).curve_potential_dam_basin(:,1)==basins_selected(i_dam).dam_height);
    i_possible = isnan(basins_selected(i_dam).curve_potential_dam_basin(:,8))==0;
    plot(basins_selected(i_dam).curve_potential_dam_basin(:,13),...
        basins_selected(i_dam).curve_potential_dam_basin(:,1),'-','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_possible,13),...
        basins_selected(i_dam).curve_potential_dam_basin(i_possible,1),'.','color',farbe(i_dam,:))
    plot(basins_selected(i_dam).curve_potential_dam_basin(i_selected,13),...
        basins_selected(i_dam).curve_potential_dam_basin(i_selected,1),'o','color',farbe(i_dam,:))
end
xlabel('share of non-suited cells (-)')
ylabel('dam height (m)')
grid on
grid minor

%% scaling and positioning of subplots

annotation('textbox',[.05,.74,.9,.23],'String','DAM CHARACTERISTICS','fontsize',schriftgroesse+2,'fontweight','bold','FitBoxToText','off')
annotation('textbox',[.05,.50,.9,.23],'String','BASIN CHARACTERISTICS','fontsize',schriftgroesse+2,'fontweight','bold','FitBoxToText','off')
annotation('textbox',[.05,.04,.9,.45],'String','DAM AND BASIN EVALUATION','fontsize',schriftgroesse+2,'fontweight','bold','FitBoxToText','off')

set(sp1,'position',[.11,.79,.22,.14],'fontsize',schriftgroesse)
set(sp2,'position',[.405,.79,.22,.14],'fontsize',schriftgroesse)
set(sp3,'position',[.7,.79,.22,.14],'fontsize',schriftgroesse)

set(sp5,'position',[.11,.55,.22,.14],'fontsize',schriftgroesse)
set(sp6,'position',[.405,.55,.22,.14],'fontsize',schriftgroesse)
set(sp7,'position',[.7,.55,.22,.14],'fontsize',schriftgroesse)

set(sp9,'position',[.11,.31,.22,.14],'fontsize',schriftgroesse)
set(sp10,'position',[.405,.31,.22,.14],'fontsize',schriftgroesse)
set(sp8,'position',[.7,.31,.22,.14],'fontsize',schriftgroesse)

set(sp11,'position',[.11,.09,.22,.14],'fontsize',schriftgroesse)
set(sp12,'position',[.405,.09,.22,.14],'fontsize',schriftgroesse)


if length(dam_ids)>1
    legend([p,p2,p3],Legende,'position',[.8,.17,.01,.01])
else
    legend([p2,p3],Legende,'position',[.8,.17,.01,.01])
end


if length(dam_ids)>1
    print('potential_dam_characteristics.pdf','-dpdf','-r300')
else
    print(sprintf('single_basin_potential_dam_characteristics_%d.pdf',i_dam),'-dpdf','-r300')
end

clear fig
function plot_5_curves_depth_storage_area(dam_ids,breite,hoehe,schriftgroesse,visible)
% ######  curves of depth, storage and area  ######
% plot curves of water depths, storage volume and flooding area to compare 
% multiple dam sites
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019

global basins_selected

clear Legende p
farbe = jet(length(basins_selected))*.8;
for i_dam = dam_ids
    Legende{i_dam} = sprintf('Basin %d',i_dam);
end

Legende{end+1} = 'dam crosssection: slope';
Legende{end+1} = 'dam crosssection: wall';


if visible == 0
    fig = figure('color',[1,1,1],'rend','painters','paperunits','centimeters','papersize',[breite,hoehe],...
        'paperpositionmode','manual','paperPosition',[0 0 breite hoehe],'position',[1 1 breite hoehe]*50,'visible','off');
else
    fig = figure('color',[1,1,1],'rend','painters','paperunits','centimeters','papersize',[breite,hoehe],...
            'paperpositionmode','manual','paperPosition',[0 0 breite hoehe],'position',[1 1 breite hoehe]*50);    
end

set(0,'CurrentFigure',fig)
subplot(1,2,1)
hold on
p1 = plot(basins_selected(dam_ids(1)).ShA_dam(:,1),basins_selected(dam_ids(1)).ShA_dam(:,2),...
    '-','color',[.5,.5,.5]);
p2 = plot(basins_selected(dam_ids(1)).ShA_wall(:,1),basins_selected(dam_ids(1)).ShA_wall(:,2),...
    '--','color',[.5,.5,.5]);
for i_dam = dam_ids
    p(i_dam) = plot(basins_selected(i_dam).ShA_dam(:,1),basins_selected(i_dam).ShA_dam(:,2),...
        '-','color',farbe(i_dam,:));
    plot(basins_selected(i_dam).ShA_wall(:,1),basins_selected(i_dam).ShA_wall(:,2),...
        '--','color',farbe(i_dam,:));
end
legend([p,p1,p2],Legende,'location','southeast','fontsize',schriftgroesse)
xlabel('basin storage volume (m³)','fontsize',schriftgroesse)
ylabel('water depth (m)','fontsize',schriftgroesse)
grid on


subplot(1,2,2)
hold on
for i_dam = dam_ids
    plot(basins_selected(i_dam).ShA_dam(:,3),basins_selected(i_dam).ShA_dam(:,2),...
        '-','color',farbe(i_dam,:));
    plot(basins_selected(i_dam).ShA_wall(:,3),basins_selected(i_dam).ShA_wall(:,2),...
        '--','color',farbe(i_dam,:));
end

xlabel('basin area (m³)','fontsize',schriftgroesse)
ylabel('water depth (m)','fontsize',schriftgroesse)
grid on


print('curves_depth_storage_area.pdf','-dpdf','-r300')

clear fig
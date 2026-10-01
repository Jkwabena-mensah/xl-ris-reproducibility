function functionIEEEFinish(fh, outPath, wCm, hCm)
%FUNCTIONIEEEFINISH  Apply the shared typography to every axis in a figure and
%export it as a vector PDF at its true printed size.
%
%Call ONCE, after all data and the legend are in place. Anything the caller has
%set deliberately (tick values, limits, legend location) is left alone; this
%only normalises type, line weights on the axes themselves, and the grid.
%
%  functionIEEEFinish(gcf, 'figures/foo.pdf')            % 8.8 x 6.6 cm
%  functionIEEEFinish(gcf, 'figures/foo.pdf', 8.8, 9.6)  % taller, e.g. 2 panels
%
%Version 1.0 (2026-09-18). License: GPLv2.

S = functionIEEEPalette();
if nargin < 3 || isempty(wCm), wCm = S.wCol; end
if nargin < 4 || isempty(hCm), hCm = 6.6;    end

set(fh,'Color','w','Units','centimeters');
pp = get(fh,'Position');
set(fh,'Position',[pp(1) pp(2) wCm hCm]);

for ax = findall(fh,'Type','axes')'
    set(ax, 'FontName',S.fnt, 'FontSize',S.fsTick, ...
            'LineWidth',S.lwAx, 'Box','on', 'TickDir','out', ...
            'XColor',[0.25 0.25 0.25], 'YColor',[0.25 0.25 0.25]);
    try, set(ax,'TickLabelInterpreter','latex'); catch, end
    grid(ax,'on');
    set(ax,'GridLineStyle',':','GridAlpha',0.22,'GridColor',[0.3 0.3 0.3]);
    %minor grids on a log axis produce a dense dotted haze that competes with
    %the data; the major decade lines carry the reading on their own
    set(ax,'XMinorGrid','off','YMinorGrid','off', ...
           'XMinorTick','off','YMinorTick','off');
    %axis labels sit at 9 pt regardless of how they were created
    set([ax.XLabel ax.YLabel],'FontName',S.fnt,'FontSize',S.fsLab, ...
        'Color',[0.10 0.10 0.10]);
end

for lg = findall(fh,'Type','legend')'
    set(lg,'FontName',S.fnt,'FontSize',S.fsLeg,'AutoUpdate','off', ...
           'Color','w','EdgeColor',[0.78 0.78 0.78],'LineWidth',0.5);
end

d = fileparts(outPath);
if ~isempty(d) && ~isfolder(d), mkdir(d); end
exportgraphics(fh, outPath, 'ContentType','vector', ...
               'BackgroundColor','white', 'Padding','figure');
fprintf('  wrote %s\n', outPath);
end

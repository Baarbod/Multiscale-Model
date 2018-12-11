function visualize(G,variable,xyz,fignum,title,clabel,climits,d)

%% ------- Figure Properties ------- %%
figure(fignum);
fig = gcf;
fig.Color = 'w';

plot(G,'XData',xyz(:,1),'YData',xyz(:,2),'ZData',xyz(:,3),...
    'EdgeCData',variable,'LineWidth',d,'EdgeAlpha',1,...
    'MarkerSize',0.1,'NodeLabel',[]);


%% ------- Axis Properties ------- %%
ax = gca;
ax.FontWeight = 'bold';
ax.FontName = 'Times New Roman';
ax.Title.String = title;
ax.Title.FontSize = 14;
ax.AmbientLightColor = 'magenta';
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.ZGrid = 'off';
ax.XAxis.Label.String = 'x (\mum)';
ax.YAxis.Label.String = 'y (\mum)';
ax.ZAxis.Label.String = 'z (\mum)';
ax.Color = 'k';

%% ------- Colormap Properties ------- %%
c = colorbar;
c.Limits = climits;
c.Label.String = clabel;
c.Location = 'southoutside';
c.FontSize = 12;
c.FontWeight = 'bold';
c.FontName = 'Times New Roman';

map = autumn;
colormap(flipud(map)) % reverse the colormap
axis equal

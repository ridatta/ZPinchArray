function formatPlots()
% Default plot properties
fnt = 'Times New Roman';
fnt_size = 24;
axis_linewidth = 2.5; 
linewidth = 2;
marker_size = 12; 


set(gcf,'Position',[0 0 800 800]); 
%set(gcf, 'Position', [pos(1) pos(2) width*100 height*100]); %<- Set size 
ax = gca;
        set(gca, 'Fontname',fnt);
        set(gca, 'Fontsize',fnt_size);
        set(gca, 'LineWidth',axis_linewidth);
        set(gca,'XMinorTick','on','YMinorTick','on')
        legend('Interpreter','latex','Fontsize',fnt_size-3);
        legend boxoff
        box on
        % axis square; 
     
end
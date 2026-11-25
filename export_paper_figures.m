function export_paper_figures(out,coarse_file,fine_file,example)
% Export the six panels in row order: I1 on top and I2 on the bottom.
% Columns: coarse reconstruction, fine reconstruction, fine signed error.
coarse=load(fullfile(out,coarse_file));
fine=load(fullfile(out,fine_file));
values={coarse.z1,fine.z1,fine.z1-fine.zz; ...
        coarse.z2,fine.z2,fine.z2-fine.zz};
suffix={'L30_dk0p5','L60_dk0p1','error_L60_dk0p1'};
folder=fullfile(out,'figures');
if ~isfolder(folder), mkdir(folder); end

if strcmp(example,'bear')
    palette=parula(256);
    low=[-.282 -.112 -.545;-.244 -.112 -.459];
    high=[.828 .910 .424;.846 .910 .435];
else
    palette=jet(256);
    recon=[coarse.z1;coarse.z2;fine.z1;fine.z2];
    low=[min(recon)*ones(2,2),[min(values{1,3});min(values{2,3})]];
    high=[max(recon)*ones(2,2),[max(values{1,3});max(values{2,3})]];
end

for i=1:2
    for j=1:3
        fig=figure('Visible','off','Color','white','Units','inches','Position',[1 1 5.2 4.6]);
        ax=axes(fig);
        imagesc(ax,fine.X(1,:),fine.Y(:,1),reshape(values{i,j},size(fine.X)));
        axis(ax,'xy'); axis(ax,'image');
        xlim(ax,[-2 2]); ylim(ax,[-2 2]);
        xticks(ax,-2:.5:2); yticks(ax,-2:.5:2);
        xlabel(ax,'X'); ylabel(ax,'Y');
        colormap(ax,palette);
        clim(ax,[low(i,j),high(i,j)]);
        colorbar(ax);
        name=sprintf('IS%d_%s.png',i,suffix{j});
        exportgraphics(fig,fullfile(folder,name),'Resolution',300,'BackgroundColor','white');
        close(fig);
    end
end
end

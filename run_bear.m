function run_bear(window_kind)
% Bear source: the two Fig. 5 settings, with and without measurement noise.
% run_bear('lanczos') repeats the same measurements with frequency tapering.
if nargin<1, window_kind='sharp'; end
window_kind=validatestring(window_kind, ...
    {'sharp','lanczos','lanczos2','hann','hann_both','cosine50'});
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'results','bear',window_kind);
if ~isfolder(out), mkdir(out); end
para=parameters();
para.window_kind=window_kind;
data=for_save_data('bear');
[A_noisy,B_noisy]=add_noise(data.U,data.B,para.noise,para.seed);

xx=-2:para.dx_inverse:2;
[X,Y]=meshgrid(xx,xx);
points=[X(:),Y(:)];
zz=bear_head(points);
edges=bear_boundary_masks(reshape(zz,size(X)),para.dx_inverse);
edge_band=edges.all{2}(:);  % Within 0.05 of any jump interface
outside=zz==0 & ~edge_band;
cases=[30 .5 30 .5;60 .1 50 .1];
rows=zeros(8,10);
row=0;
for noise=[0,para.noise]
    for c=1:size(cases,1)
        L=cases(c,1); dk=cases(c,4);
        k=cases(c,2):dk:cases(c,3);
        sensors=1:(numel(data.theta)/L):numel(data.theta);
        frequencies=round(k/data.k_total(1));
        theta=data.theta(sensors);
        if noise==0
            A=data.U(sensors,frequencies); B=data.B(sensors,frequencies);
        else
            A=A_noisy(sensors,frequencies); B=B_noisy(sensors,frequencies);
        end
        para.lambda=dk:dk:para.lambda_max;
        [z1,z2,check]=Simple_f_yy(points,A,B,k,theta,para);
        assert(all(isfinite(z1)) && all(isfinite(z2)));
        assert(max(check.inverse1_interp_max,check.inverse2_interp_max)<1e-7);

        error_map=[z1-zz,z2-zz];
        for j=1:2
            err=error_map(:,j);
            row=row+1;
            rows(row,:)=[cases(c,:),noise,j,norm(err)/norm(zz), ...
                norm(err(edge_band))/norm(zz(edge_band)), ...
                sqrt(mean(err(edge_band).^2)),sqrt(mean(err(outside).^2))];
        end
        file=sprintf('case_%d_noise%d.mat',c,round(100*noise));
        save(fullfile(out,file),'X','Y','zz','z1','z2','k','theta','para','noise','check');
    end
end
T=array2table(rows,'VariableNames', ...
    {'L','k_min','k_max','dk','noise','indicator','RE','near_jump_RE','near_jump_RMSE','outside_RMS'});
writetable(T,fullfile(out,'metrics.csv'));
disp(T);
export_paper_figures(out,sprintf('case_1_noise%d.mat',round(100*para.noise)), ...
    sprintf('case_2_noise%d.mat',round(100*para.noise)),'bear');
end

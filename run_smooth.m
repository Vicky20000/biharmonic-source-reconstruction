function run_smooth()
% Smooth source: four Table 1 cases and two noise-free reconstructions.
root = fileparts(mfilename('fullpath'));
out = fullfile(root,'results','smooth');
if ~isfolder(out), mkdir(out); end
para = parameters();
data = for_save_data('smooth');
[A_noisy,B_noisy] = add_noise(data.U,data.B,para.noise,para.seed);

xx = -2:para.dx_inverse:2;
[X,Y] = meshgrid(xx,xx);
points = [X(:),Y(:)];
zz = general_source_function(points);

% Columns: number of sensors, k_min, k_max, dk, relative noise level.
cases = [30 .5 30 .5 para.noise;
         60 .5 30 .5 para.noise;
         30 .1 50 .1 para.noise;
         60 .1 50 .1 para.noise;
         30 .5 30 .5 0;
         60 .1 50 .1 0];
RE = zeros(size(cases,1),2);
for c=1:size(cases,1)
    L=cases(c,1); dk=cases(c,4); noise=cases(c,5);
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
    RE(c,:)=[norm(z1-zz),norm(z2-zz)]/norm(zz);
    save(fullfile(out,sprintf('case_%02d.mat',c)), ...
        'X','Y','zz','z1','z2','k','theta','para','noise','check');
end

T=array2table([cases,RE],'VariableNames', ...
    {'L','k_min','k_max','dk','noise','RE_I1','RE_I2'});
writetable(T,fullfile(out,'summary.csv'));
Indicator={'I_S^(1)';'I_S^(2)';'I_S^(1)';'I_S^(2)'};
L=[30;30;60;60];
RE_dk05=[RE(1,:).';RE(2,:).'];
RE_dk01=[RE(3,:).';RE(4,:).'];
Table1=table(Indicator,L,RE_dk05,RE_dk01);
writetable(Table1,fullfile(out,'table1_relative_errors.csv'));
disp(Table1);
export_paper_figures(out,'case_01.mat','case_04.mat','smooth');
end

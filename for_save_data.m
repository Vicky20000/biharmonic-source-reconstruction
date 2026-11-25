function data = for_save_data(example,rebuild)
% Save the scattered field and its Laplacian on a common fine sensor grid.
% Set rebuild=true after changing the source function or forward kernel.
if nargin<2, rebuild=false; end
example = validatestring(example,{'smooth','bear'});
root = fileparts(mfilename('fullpath'));
folder = fullfile(root,'data');
if ~isfolder(folder), mkdir(folder); end
file = fullfile(folder,[example '_forward.mat']);

para = parameters();
theta = (0:59)*2*pi/60;
k_total = 0.1:0.1:50;
z = para.radius*[cos(theta(:)),sin(theta(:))];
forward_settings = [para.dx_forward,para.distance_step];
if isfile(file) && ~rebuild
    saved = load(file);
    if isequal(saved.z,z) && isequal(saved.k_total,k_total) && ...
            isequal(saved.forward_settings,forward_settings)
        data = saved;
        fprintf('Loaded %s\n',file);
        return
    end
end

source = @general_source_function;
if strcmp(example,'bear'), source=@bear_head; end
[U,B,info]=u_s(z,k_total,para.dx_forward,para.distance_step,source);
data=struct('U',U,'B',B,'theta',theta,'z',z,'k_total',k_total, ...
    'forward_settings',forward_settings,'info',info);
save(file,'-struct','data');
fprintf('Saved %s\n',file);
end

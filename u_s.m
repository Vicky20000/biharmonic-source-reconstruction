function [U_s_R2, Laplace_U_s_R2, info] = u_s(z_matrix,k_matrix,dx,distance_step,source)
% Evaluate u^s and Delta u^s by composite midpoint quadrature on [-3,3]^2.
% Rows correspond to sensors and columns to wavenumbers. Set distance_step
% to zero to evaluate the kernels directly at the quadrature distances.
if nargin<3, dx=0.01; end
if nargin<4, distance_step=0.001; end
if nargin<5, source=@general_source_function; end
x_range = -3+dx/2:dx:3-dx/2;
[X,Y] = meshgrid(x_range,x_range);
y = [X(:),Y(:)];
source_values = source(y);
% Only midpoints in the support contribute.
y = y(source_values~=0,:);
source_values = source_values(source_values~=0)*dx^2;
num_z = size(z_matrix,1); num_k = numel(k_matrix);
U_s_R2 = complex(zeros(num_z,num_k));
Laplace_U_s_R2 = U_s_R2;
if distance_step>0
    % Distribute the midpoint weights to four neighboring distance nodes.
    % This is cubic Lagrange interpolation of the radial Bessel kernels.
    d = (0:distance_step:7).';
    W = zeros(numel(d),num_z);
    for j=1:num_z
        dist = hypot(y(:,1)-z_matrix(j,1),y(:,2)-z_matrix(j,2));
        t = dist/distance_step; n=floor(t); a=t-n;
        weights = [-a.*(a-1).*(a-2)/6, (a+1).*(a-1).*(a-2)/2,...
                   -(a+1).*a.*(a-2)/2, (a+1).*a.*(a-1)/6];
        for q=1:4
            idx = n+q-1; % MATLAB index for distance nodes n-1 through n+2.
            W(:,j)=W(:,j)+accumarray(idx,source_values.*weights(:,q),[numel(d),1]);
        end
    end
    keep=any(W~=0,2); d=d(keep); W=W(keep,:);
end
for m=1:num_k
    k=k_matrix(m);
    if distance_step>0
        H=besselh(0,1,k*d); K=besselk(0,k*d);
        % Fundamental solution and its Laplacian at the sensor.
        U_s_R2(:,m)=W.'*(1i/(8*k^2)*(H+2i/pi*K));
        Laplace_U_s_R2(:,m)=W.'*(-1/8*(1i*H+2/pi*K));
    else
        for j=1:num_z
            dist=hypot(y(:,1)-z_matrix(j,1),y(:,2)-z_matrix(j,2));
            H=besselh(0,1,k*dist); K=besselk(0,k*dist);
            U_s_R2(j,m)=sum(source_values.*(1i/(8*k^2)*(H+2i/pi*K)));
            Laplace_U_s_R2(j,m)=sum(source_values.*(-1/8*(1i*H+2/pi*K)));
        end
    end
end
info.dx=dx; info.full_grid_size=[numel(x_range),numel(x_range)];
info.nonzero_nodes=numel(source_values); info.distance_step=distance_step;
info.quadrature='composite midpoint';
end

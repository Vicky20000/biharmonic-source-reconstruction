function [f_yy_1,f_yy_2,check]=Simple_f_yy(yy,A,B,k,theta,para)
% Reconstruct both indicators at the vertices yy using the fields A and B.
% k, r and lambda use left rectangles. The L angles have weight 2*pi/L.
R=para.radius; L=numel(theta); wk=left_weights(k);
if ~isfield(para,'window_kind'), para.window_kind='sharp'; end
window_k=spectral_window(k,k(end),para.window_kind);
wk2=wk.*window_k; % Taper the measurement frequency k for I2.
r=para.r; wr=left_weights(r);
% Use the same step for k and lambda, with a separate cutoff for lambda.
if isempty(para.lambda)
    dlambda=k(2)-k(1); lambda=dlambda:dlambda:para.lambda_max;
else
    lambda=para.lambda;
end
% Integrate from the first positive lambda node to the cutoff.
window_lambda=spectral_window(lambda,lambda(end),para.window_kind);
wlambda=left_weights(lambda).*window_lambda; % Taper the outer lambda integral for I1.
% For I1, sum first over k, then over r. Each row is one sensor.
inner=(imag(A).*(wk.*k.^3))*besselj(0,k(:)*r);
C=(inner.*(wr.*r))*besselj(0,r(:)*lambda);
D=(inner.*(wr.*r))*bessely(0,r(:)*lambda);
C=C.*(wlambda.*lambda.^2); D=D.*(wlambda.*lambda.^2);
% Real part of I2, written using the J1 and Y1 kernels.
E=(k.^2.*real(A)-real(B)).*(wk2.*k.^2);
F=imag(A).*(2*wk2.*k.^4);
d=(0.16:para.distance_step:5.84).';
J1=besselj(1,d*lambda); Y1=bessely(1,d*lambda);
T1=Y1*C.'-J1*D.';
J1=besselj(1,d*k); Y1=bessely(1,d*k);
T2=J1*E.'+Y1*F.';
f_yy_1=zeros(size(yy,1),1); f_yy_2=f_yy_1;
% Evaluate the radial kernels directly at 41 points to check interpolation.
idx=unique(round(linspace(1,size(yy,1),41)));
ref1=zeros(numel(idx),1); ref2=ref1;
for j=1:L
    x=R*[cos(theta(j)),sin(theta(j))];
    delta=yy-x;
    dist=hypot(delta(:,1),delta(:,2));
    factor=(delta(:,1)*cos(theta(j))+delta(:,2)*sin(theta(j)))./dist;
    f_yy_1=f_yy_1+factor.*interp1(d,T1(:,j),dist,'spline');
    f_yy_2=f_yy_2+factor.*interp1(d,T2(:,j),dist,'spline');
    dd=dist(idx);
    ref1=ref1+factor(idx).*(bessely(1,dd*lambda)*C(j,:).'-besselj(1,dd*lambda)*D(j,:).');
    ref2=ref2+factor(idx).*(besselj(1,dd*k)*E(j,:).'+bessely(1,dd*k)*F(j,:).');
end
f_yy_1=f_yy_1*(R/pi)*(2*pi/L);
f_yy_2=f_yy_2*(R/(2*pi))*(2*pi/L);
ref1=ref1*(R/pi)*(2*pi/L); ref2=ref2*(R/(2*pi))*(2*pi/L);
check.quadrature='left endpoint rectangle';
check.window_k=window_k; check.window_lambda=window_lambda;
check.window_kind=para.window_kind;
check.lambda_nodes=lambda;
check.k_upper_weight=wk(end); check.r_upper_weight=wr(end); check.lambda_upper_weight=wlambda(end);
check.inverse1_interp_max=max(abs(ref1-f_yy_1(idx)));
check.inverse2_interp_max=max(abs(ref2-f_yy_2(idx)));
end

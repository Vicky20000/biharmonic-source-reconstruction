function w=spectral_window(q,K,kind)
% Weights on the positive-frequency half of a symmetric window.
% The Lanczos window uses sinc(q/K), with sinc(0)=1.
% cosine50 has a flat passband up to K/2.
t=q/K; w=ones(size(q));
switch kind
 case 'sharp'
 case {'lanczos','lanczos2'}
    use=abs(t)>eps; w(use)=sin(pi*t(use))./(pi*t(use));
    if strcmp(kind,'lanczos2'), w=w.^2; end
 case {'hann','hann_both'}
    w=.5*(1+cos(pi*t));
    if strcmp(kind,'hann_both')
        dq=q(2)-q(1); u=max(0,min(1,(q-q(1))/(2*dq)));
        w=w.*(.5-.5*cos(pi*u)); % Taper the first two intervals at the low end.
    end
 case 'cosine50'
    u=max(0,min(1,(q-.5*K)/(.5*K))); w=.5*(1+cos(pi*u));
 otherwise, error('Unknown window.');
end
w(q>K)=0;
assert(all(isfinite(w)) && all(w>=-1e-14 & w<=1+1e-14),'Invalid window weights.');
end

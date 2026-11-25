function [A,B,noise_matrix]=add_noise(U,V,noise,seed)
% Multiply each real and imaginary component by 1 + noise*xi.
% Use four independent arrays of uniform draws on [-1,1].
rng(seed,'twister');
noise_matrix=2*rand([size(U),4])-1;
A=real(U).*(1+noise*noise_matrix(:,:,1))+1i*imag(U).*(1+noise*noise_matrix(:,:,2));
B=real(V).*(1+noise*noise_matrix(:,:,3))+1i*imag(V).*(1+noise*noise_matrix(:,:,4));
end

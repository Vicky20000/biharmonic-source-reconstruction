function para = parameters()
% Common settings for the smooth and piecewise-constant source examples.
para.radius = 3;
para.dx_forward = 0.01;    % Midpoint grid on [-3,3]^2
para.dx_inverse = 0.01;    % Vertex grid on [-2,2]^2
para.noise = 0.20;
para.seed = 20261120;
para.lambda_max = 30;
para.lambda = [];         % Use the same step as k by default.
para.dr = 0.01;
para.r = para.dr:para.dr:6;
para.distance_step = 0.001;
para.window_kind = 'sharp';
end

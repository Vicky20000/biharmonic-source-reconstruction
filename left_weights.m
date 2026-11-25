function w=left_weights(x)
% Left rectangle weights on [x(1), x(end)], including nonuniform grids.
% x(end) is the end of the interval and has zero weight.
x=x(:).';
assert(numel(x)>=2 && all(isfinite(x)) && all(diff(x)>0),'Quadrature nodes must be finite and strictly increasing.');
w=[diff(x),0];
end

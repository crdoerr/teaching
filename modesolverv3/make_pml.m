function s = make_pml(coord, pml_thick, sigma_max, order, k0)
cmin = min(coord);
cmax = max(coord);
s = ones(size(coord));

left = coord < cmin + pml_thick;
right = coord > cmax - pml_thick;

dl = (cmin + pml_thick - coord(left))/pml_thick;
dr = (coord(right) - (cmax - pml_thick))/pml_thick;

sigma = zeros(size(coord));
sigma(left) = sigma_max * dl.^order; % dampening profile for left bound.
sigma(right) = sigma_max * dr.^order; % dampening profile for right bound.

s = 1 - 1i*sigma/k0; % complex coordinate stretching for boundary (no reflections)
end

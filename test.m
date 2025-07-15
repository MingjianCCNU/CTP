n = 6;
m = 5;
c    = cell(1, n);
c(:) = {1:m};
[d{1:n}] = ndgrid(c{:});
allp = reshape(cat(n+1, d{:}), [], n);
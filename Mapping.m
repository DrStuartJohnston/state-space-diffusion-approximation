%% Function to map Cartesian co-ordinates to a vector index

function index = Mapping(dimensions,indices)

nDim = numel(dimensions);
index = 0;
for i = nDim:-1:2
    index = index + (indices(i)-1)*prod(dimensions(1:i-1));
end
index = index + indices(1);
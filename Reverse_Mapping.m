%% Function to map a vector index to Cartesian co-ordinates

function coord = Reverse_Mapping(dimensions,index)

nDim = numel(dimensions);
coord = zeros(nDim,1);
for i = nDim:-1:1
   coord(i) = ceil(index/prod(dimensions(1:i-1)));
   index = index - (coord(i)-1)*prod(dimensions(1:i-1));
end
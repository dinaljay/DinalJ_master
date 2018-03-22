function vector = vector(matrix)
%vector Makes the matrix 1xn vector matrix
%   Input = nD matrix

vector = reshape(matrix,[1,numel(matrix)]);
end


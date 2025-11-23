function [bestIndices, minNorm] = brute_force(A, k)
% BRUTE_FORCE  Selects the k-column subset of A whose pseudoinverse 
%              has the smallest Frobenius norm.
%
% INPUTS:
%   A : m-by-n matrix
%   k : number of columns to choose
%
% OUTPUTS:
%   bestIndices : the indices of the selected k columns
%   minNorm     : the minimal Frobenius norm of the pseudoinverse among all subsets

    % Get matrix dimensions
    [m, n] = size(A);

    % Validate that k is not larger than the number of columns
    if k > n
        error('k must be less than or equal to the number of columns in A');
    end

    % Generate all possible combinations of k columns from 1:n
    combos = nchoosek(1:n, k);
    numCombos = size(combos, 1);

    % Initialize tracking variables:
    %   minNorm     = best (smallest) Frobenius norm found so far
    %   bestIndices = column indices corresponding to that norm
    minNorm = inf;
    bestIndices = [];

    % Loop over all possible combinations
    for i = 1:numCombos
        % Get the i-th set of column indices
        cols = combos(i, :);

        % Extract the submatrix consisting of those k columns
        submatrix = A(:, cols);

        % Compute its pseudoinverse
        pinvSub = pinv(submatrix);

        % Compute the Frobenius norm of the pseudoinverse
        normVal = norm(pinvSub, 'fro');

        % If this subset yields a smaller norm, update the best-so-far
        if normVal < minNorm
            minNorm = normVal;
            bestIndices = cols;
        end
    end
end

function [selected_matrix, row_indices] = luke_algorithm(A, k)
% LUKE_ALGORITHM   Selects k rows of A using the LUKE greedy algorithm.
%
%   The algorithm iteratively removes the row whose removal
%   maximizes the smallest singular value of the remaining matrix.
%
% INPUTS:
%   A : n-by-m matrix
%   k : number of rows to select
%
% OUTPUTS:
%   selected_matrix : the resulting k-by-m matrix after selection
%   row_indices     : the indices of the selected rows

    % Number of rows in the matrix
    [n, ~] = size(A);

    % Validate that k is not larger than the number of available rows
    if k > n
        error('k must be less than or equal to the number of rows in the matrix.');
    end

    % Start with all rows selected
    row_indices = 1:n; 
    
    % Continue removing rows until only k remain
    while length(row_indices) > k

        % Track the best (largest) minimal singular value after a removal
        best_sigma_min = -inf;

        % Track which row should be removed
        worst_row_to_remove = -1;
        
        % Try removing each row one at a time
        for i = 1:length(row_indices)
         
            % Remove row i hypothetically
            temp_indices = row_indices;
            temp_indices(i) = [];
            
            % Build the temporary submatrix without that row
            temp_matrix = A(temp_indices, :);
     
            % Compute the smallest singular value of the remaining matrix
            sigma_min = min(svd(temp_matrix));
            
            % The LUKE rule:
            % Remove the row whose removal yields the *largest* σ_min
            if sigma_min > best_sigma_min
                best_sigma_min = sigma_min;
                worst_row_to_remove = i;
            end
        end
        
        % Permanently remove the worst row
        row_indices(worst_row_to_remove) = [];
    end
    
    % Final selected submatrix
    selected_matrix = A(row_indices, :);
end

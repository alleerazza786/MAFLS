function computeAndDisplayResultsForDifferentN()
    % Fixed value for m
    m = 2;
    % Varying values for n from 2 to 15
    n_values = 2:15;
    % Initialize arrays to store the results
    n_array = [];    SRAE_array = [];
    PAE_array = [];    APL_array = [];
    NTS_array = [];    GMPT_array = [];
    KRI_array = [];    GE_array = [];   SR_array = [];
    for n = n_values
        % Call the function to compute results
        [SRAE, PAE, APL, NTS, GMPT, KRI] = computeResults(m, n);        
        % Call the function to evaluate GE
        GE_result = evaluate_GE(m, n);        
        % Call the function to calculate SR
        SR_result = calculate_SR(m, n);
        % Add the results to the arrays
        n_array = [n_array; n];
        SRAE_array = [SRAE_array; SRAE];
        PAE_array = [PAE_array; PAE];
        APL_array = [APL_array; APL];
        NTS_array = [NTS_array; NTS];
        GMPT_array = [GMPT_array; GMPT];
        KRI_array = [KRI_array; KRI];
        GE_array = [GE_array; GE_result];
        SR_array = [SR_array; SR_result];
    end
    % Create a table with the results
    resultsTable = table(n_array, APL_array,  GMPT_array, KRI_array, GE_array, SR_array, ...
        'VariableNames', {'n', 'APL', 'GMPT', 'KRI', 'GE', 'SR'});
    % Display the table
    disp(resultsTable);
end
function [SRAE, PAE, APL, NTS, GMPT, KRI] = computeResults(m, n)
    % Check if the inputs are positive integers
    if ~isscalar(m) || ~isnumeric(m) || m <= 0 || mod(m, 1) ~= 0 || ...
            ~isscalar(n) || ~isnumeric(n) || n <= 0 || mod(n, 1) ~= 0
        error('Inputs must be positive integers.');
    end
    % Calculate the adjacency eigenvalues using the provided formula
    adjacencyEigenvalues = zeros(m * n, 1);
    laplacianEigenvalues = zeros(m * n, 1);
    signlessLaplacianEigenvalues = zeros(m * n, 1);
    idx = 1;
    for lambda = 0:(m-1)
        for nu = 0:(n-1)
            adjacencyEigenvalues(idx) = 2 * cos(pi * lambda / (m+1)) + 2 * cos(pi * nu / (n+1));
            laplacianEigenvalues(idx) = 4 - 2 * cos(pi * lambda / m) - 2 * cos(pi * nu / n);
            signlessLaplacianEigenvalues(idx) = 4 + 2 * cos(2 * pi * lambda / (m + 1)) + 2 * cos(2 * pi * nu / (n + 1));
            idx = idx + 1;
        end
    end
    % Calculate the SRAE (Sum of Reciprocals of All Laplacian Eigenvalues, excluding zero)
    SRAE = sum(1./laplacianEigenvalues(laplacianEigenvalues ~= 0));
    % Calculate the PAE (Product of All Laplacian Eigenvalues, excluding zero)
    PAE = prod(laplacianEigenvalues(laplacianEigenvalues ~= 0));
    % Calculate the APL (Average Path Length)
    APL = SRAE * m * n / (2 * (m * n - 1));
    % Calculate the NTS
    NTS = PAE * 2 / (m * n);
    % Calculate the GMPT
    GMPT = ((2 * (2 * m * n - m - n)) / (m * n * (m * n - 1))) * SRAE;
    % Calculate the KRI
    KRI = m * n * SRAE;
end
function GE_value = evaluate_GE(m, n)
    % Initialize the sum
    GE_value = 0;    
    % Nested loops to evaluate the double summation
    for lambda = 0:(m-1)
        for nu = 0:(n-1)
            % Calculate the inner expression
            inner_expression = abs(4 - 2*cos(pi*lambda/m) - 2*cos(pi*nu/n) - (2*(2*m*n - m - n))/(m*n));            
            % Add the inner expression to the sum
            GE_value = GE_value + inner_expression;
        end
    end
end
function SR_value = calculate_SR(m, n)
    % Initialize the max value
    max_value = -Inf;    
    % Nested loops to find the maximum value of the expression
    for lambda = 0:(m-1)
        for nu = 0:(n-1)
            % Calculate the expression
            expression_value = abs(4 - 2*cos(pi*lambda/m) - 2*cos(pi*nu/n));        
            % Update max_value if needed
            max_value = max(max_value, expression_value);
        end
    end   
    % Store the max value for the current n
    SR_value = max_value;
end

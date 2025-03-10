% Lab 4.1 --> Big M method
clc; clear; close all;

% Input basic details
no_of_var = input("Enter the number of variables: ");
no_of_constraints = input("Enter the number of constraints: ");
LessThanEqualTo = input("Enter the no of less than equal to constraints: ");
EqualTo = input("Enter the no of equal to constraints: ");
GreaterThanEqualTo = input("Enter the no of greater than equal to constraints: ");

% Input matrices
A = input("Enter the matrix A: ");
b = input("Enter the constant matrix (RHS): ");
c = input("Enter the coefficients of the objective function: ");

% Define a large M value
M = 1e6;

% Initialize extra columns for slack, surplus, and artificial variables
extra_mat = zeros(no_of_constraints, LessThanEqualTo + EqualTo + 2 * GreaterThanEqualTo);
extra_mat_for_objective = zeros(1, size(extra_mat, 2));
CC = [c, zeros(1,no_of_constraints-EqualTo)];
% Add slack variables (for <= constraints)
for i = 1:LessThanEqualTo
    extra_mat(i, i) = 1;
end

% Add artificial variables (for = constraints)
for i = 1:EqualTo
    extra_mat(LessThanEqualTo + i, LessThanEqualTo + i) = 1;
    extra_mat_for_objective(LessThanEqualTo + i) = -1;
end

% Add surplus and artificial variables (for >= constraints)
for i = 1:GreaterThanEqualTo
    index = LessThanEqualTo + EqualTo + i;
    extra_mat(index, index) = -1;
    extra_mat(index, index + 1) = 1;
    extra_mat_for_objective(index + 1) = -1;
end

% Update A and c with extra variables
A = [A, extra_mat];
C = [zeros(size(c)), extra_mat_for_objective];
disp("Objective Function Coefficients in Phase 1:");
disp(C);

% Initialize Simplex Table
table = [A, b];
Cb = zeros(1, no_of_constraints);
Cb(LessThanEqualTo + 1:end) = -1;

Zj_minus_Cj = Cb * A - C;
disp("Initial Simplex Table:");
disp(table);

% Simplex Iteration
no_of_iter = 0;
max_iterations = nchoosek(no_of_var + no_of_constraints, no_of_constraints);

while any(Zj_minus_Cj < 0)
    no_of_iter = no_of_iter + 1;
    if no_of_iter > max_iterations
        disp("No solution exists.");
        return;
    end
    
    % Find entering variable
    [~, entering_col] = min(Zj_minus_Cj);
    
    % Compute ratios for pivot row
    ratios = inf(no_of_constraints, 1);
    valid_rows = table(:, entering_col) > 0;
    ratios(valid_rows) = table(valid_rows, end) ./ table(valid_rows, entering_col);
    
    % Check for unbounded solution
    if all(ratios == inf)
        disp("The solution is unbounded.");
        return;
    end
    
    % Find pivot row
    [~, pivot_row] = min(ratios);
    fprintf("Pivot element at row %d, column %d.\n", pivot_row, entering_col);
    
    % Normalize pivot row
    table(pivot_row, :) = table(pivot_row, :) / table(pivot_row, entering_col);
    
    % Update other rows
    for i = 1:no_of_constraints
        if i ~= pivot_row
            table(i, :) = table(i, :) - table(i, entering_col) * table(pivot_row, :);
        end
    end
    
    % Update Cb values
    Cb(pivot_row) = C(entering_col);
    
    % Recalculate Zj - Cj
    Zj_minus_Cj = Cb * table(:, 1:end-1) - C;
    
    % Display updated table
    fprintf("\nIteration %d:\n", no_of_iter);
    disp("Updated Simplex Table:");
    disp(table);
    disp("Zj - Cj:");
    disp(Zj_minus_Cj);
end

% Extract optimal solution
disp("Optimal Solution:");
X = zeros(1, no_of_var + no_of_constraints);
for i = 1:no_of_constraints
    basic_var_index = find(table(i, 1:no_of_var + no_of_constraints) == 1);
    if ~isempty(basic_var_index)
        X(basic_var_index) = table(i, end);
        if(C(1,basic_var_index)==(-1))
            if(X(basic_var_index)~=0)
                disp("NO solution");
               return;
            end
        end
    end
end
% Display results
disp("Decision Variables:");
disp(X);
disp("Optimal Objective Value:");
disp(Cb * table(:, end));
[m,n] = size(C);
CCB = zeros(0);
table1 = zeros(0);
for i  = 1:n
    if(C(1,i)==0)
        for j = 1:no_of_constraints
            if(table(j,i)==1)
                CCB(j,1) = CC(1,i);
            end
        end
        table1 = [table1,table(1:no_of_constraints,i)];
    end
end
table1 = [table1 table(1:no_of_constraints,end)];
CCB = CCB';
 Zj_minus_Cj = CCB * table1(:, 1:end-1) - CC;

disp("Phase 2");
while any(Zj_minus_Cj < 0)
    no_of_iter = no_of_iter + 1;
    if no_of_iter > max_iterations
        disp("No solution exists.");
        return;
    end
    
    % Find entering variable
    [~, entering_col] = min(Zj_minus_Cj);
    
    % Compute ratios for pivot row
    ratios = inf(no_of_constraints, 1);
    valid_rows = table1(:, entering_col) > 0;
    ratios(valid_rows) = table(valid_rows, end) ./ table(valid_rows, entering_col);
    
    % Check for unbounded solution
    if all(ratios == inf)
        disp("The solution is unbounded.");
        return;
    end
    
    % Find pivot row
    [~, pivot_row] = min(ratios);
    fprintf("Pivot element at row %d, column %d.\n", pivot_row, entering_col);
    
    % Normalize pivot row
    table1(pivot_row, :) = table1(pivot_row, :) / table1(pivot_row, entering_col);
    
    % Update other rows
    for i = 1:no_of_constraints
        if i ~= pivot_row
            table1(i, :) = table1(i, :) - table1(i, entering_col) * table1(pivot_row, :);
        end
    end
    
    % Update Cb values
    CCB(pivot_row) = CC(entering_col);
    
    % Recalculate Zj - Cj
    Zj_minus_Cj = CCB * table1(:, 1:end-1) - CC;
    
    % Display updated table
    fprintf("\nIteration %d:\n", no_of_iter);
    disp("Updated Simplex Table:");
    disp(table1);
    disp("Zj - Cj:");
    disp(Zj_minus_Cj);
end
disp("Optimal Solution:");
disp("Decision Variables:");
disp(X);
disp("Optimal Objective Value:");
disp(CCB * table1(:, end));
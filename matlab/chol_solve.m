function [rows,memory_delta,solve_time,relative_error] = chol_solve(filename, debug)
%SOLVE Summary of this function goes here
%   Detailed explanation goes here

if (nargin == 1)
    debug = 0;
end

[user] = memory;
proc_memory_start = user.MemUsedMATLAB ; % memoria iniziale

if debug
    disp('Reading input file:' + filename); 
end
load(fullfile('', 'matrix_mat', filename), "Problem");

rows = size(Problem.A, 1);

if debug
    disp('Fill solution with ones');
end
x_es = ones(size(Problem.A, 1), 1); % creo vettori 1s

if debug
    disp('Calculate b');
end
b = Problem.A*x_es;

if debug
    disp('Calculate Cholesky');
end

tic
R = chol(Problem.A);
x_ap = R\(R'\b); % cholesky.matlab official documentation
solve_time = toc;

relative_error = norm(x_ap - x_es)/norm(x_es);
if debug
    disp(['Relative error: ', num2str(relative_error)]);
end

[user] = memory;
memory_delta = user.MemUsedMATLAB - proc_memory_start;
% disp('Total mem: ' + proc_memory_end); % memoria usata
end
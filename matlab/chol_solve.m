function [rows,memory_delta,solve_time,relative_error] = chol_solve(filename, debug)
%SOLVE Summary of this function goes here
%   Detailed explanation goes here

if (nargin == 1)
    debug = 0;
end

if debug
    disp('Reading input file:' + filename);
end
load(fullfile('', 'matrix_mat', filename), "Problem");

if isunix
    [~, pid] = system('pgrep MATLAB');
    [~, mem] = system(['smem | grep ' strtrim(pid) ' | grep -v grep | awk -F "[ ]+" ''{print $6}''']);

    proc_memory_start = str2num(mem) * 1024;
else
    [user] = memory;
    proc_memory_start = user.MemUsedMATLAB;
end

rows = size(Problem.A, 1);
if debug
    disp('Fill solution with ones');
end
x_es = ones(rows, 1);
if debug
    disp('Calculate b');
end
b = Problem.A*x_es;

if debug
    disp('Calculate Cholesky');
end

tic;
R = chol(Problem.A);
x_ap = R\(R'\b); % cholesky.matlab official documentation
solve_time = toc;

if isunix
    [~, mem] = system(['smem | grep ' strtrim(pid) ' | grep -v grep | awk -F "[ ]+" ''{print $6}''']);

    proc_memory_end = str2num(mem) * 1024;
    memory_delta = proc_memory_end - proc_memory_start;
else
    [user] = memory;
    memory_delta = user.MemUsedMATLAB - proc_memory_start;
end

if debug
    disp(['Total mem: ', num2str(memory_delta)]); % memoria usata
end

relative_error = norm(x_ap - x_es)/norm(x_es);
if debug
    disp(['Relative error: ', num2str(relative_error)]);
end
end
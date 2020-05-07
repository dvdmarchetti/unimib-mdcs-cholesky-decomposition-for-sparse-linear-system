
if ~(isunix || ispc)
    disp('Platform not supported');
    exit;
end

profile clear;
profile -memory on;

path = fullfile('matrix_mat', '*.mat');
files = dir(path);
filesCount = length(files);

results = ["filename", "size", "memory_delta", "chol_time", "relative_error"];

for K = 1 : filesCount
    filename = convertCharsToStrings(files(K).name);
    [size, memory_delta, solve_time, relative_error] = chol_solve(filename, 1);
    
    results = [results(1:K,:); [strtok(filename, '.'),size,memory_delta,solve_time,relative_error]];
    clearvars s ms me bt ct re;
end

profile off;
if isunix
    writematrix(results, 'unix-output.csv');
elseif ispc
    writematrix(results, 'windows-output.csv');
end


path = fullfile('matrix', '*.mat');
files = dir(path);
for K = 1 : length(files)
    load(fullfile('matrix', files(K).name), "Problem");
    mmwrite_symmetric(fullfile('matrix_converted', strcat(files(K).name,'.mtx')), Problem.A);
end
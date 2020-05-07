path = fullfile('matrix_mat', '*.mat');
files = dir(path);

for K = 1 : length(files)
    load(fullfile('matrix_mat', files(K).name), "Problem");
    parts = strtok(files(K).name, '.');
    newname = strcat(parts,'.mtx');
    mmwrite_symmetric(fullfile('matrix_mtx', newname), Problem.A);
end
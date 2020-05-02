disp('Reading input');  % printf
load(fullfile('..', 'matrix_matlab', 'ex15.mat'), "Problem");   % carico matrice
% creo variabile Problem con il file all'interno
disp('Calculate Cholesky');      
R = chol(Problem.A);    % matrice (file senza descrizioni)

disp('Fill solution with ones');
x_es = ones(size(Problem.A, 1), 1); % creo vettori 1s

disp('Calculate b');
b = Problem.A*x_es;     % vettore termini noti

disp('Find solution');
x_ap = R\(R'\b);    % cholesky.matlab official documentation

err = norm(x_es - x_ap)/norm(x_es);
disp(['Relative error: ', num2str(err)]);   %errore relativo

% cholesky_benchmark(fullfile('..', 'matrix_matlab'));

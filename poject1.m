 % Read in necessary data
x = importdata('data_points.txt'); x = [0; x; 1];
v = importdata('numerical.txt');   v = [0; v; 0]; %Dirichlet boundary cond.
u = importdata('analytical.txt');  u = [0; u; 0]; %Dirichlet boundary cond.

 % Plot
plot(x, v, 'b-', x, u, 'r-')
title(sprintf('Number of grid points n = %d.', numel(x) - 2))
grid on
legend('numerical', 'analytical')

 % Calculate maximum relative error, and write to command window
eps = log10(abs((v - u)./u));
maxerr = max(eps);
fprintf('n = %d || Maximum relative error = %d \n', numel(x) - 2, maxerr);
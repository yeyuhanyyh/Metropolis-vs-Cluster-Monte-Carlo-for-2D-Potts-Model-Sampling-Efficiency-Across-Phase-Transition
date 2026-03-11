% Calculate the change in Hamiltonian
function deltaH = calculateDeltaH(value, sigma, i, j, J, h)
% Get the size of the lattice
[N, M] = size(sigma);

% Calculate the original Hamiltonian
H0 = -J * ((sigma(i, j) == sigma(mod(i-2, N)+1, j)) + (sigma(i, j) == sigma(mod(i, N)+1, j)) + (sigma(i, j) == sigma(i, mod(j-2, M)+1)) + (sigma(i, j) == sigma(i, mod(j, M)+1)));
H1 = -J * ((value == sigma(mod(i-2, N)+1, j)) + (value == sigma(mod(i, N)+1, j)) + (value == sigma(i, mod(j-2, M)+1)) + (value == sigma(i, mod(j, M)+1)));
H0 = H0 - h*sigma(i, j);
H1 = H1 - h*value;

% Calculate the change in Hamiltonian
deltaH = (H1-H0);
end
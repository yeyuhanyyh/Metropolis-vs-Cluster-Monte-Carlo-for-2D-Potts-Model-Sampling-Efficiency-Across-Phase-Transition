% Function to calculate the Hamiltonian
function result = calculateH(sigma, J, h)
% Get the size of the lattice
[N, M] = size(sigma);

% Initialize the Hamiltonian
H = zeros(N, M);

% Calculate the Hamiltonian
for i = 1:N
    for j = 1:M
        H(i, j) = -J * ((sigma(i, j) == sigma(mod(i-2, N)+1, j)) + (sigma(i, j) == sigma(mod(i, N)+1, j)) + (sigma(i, j) == sigma(i, mod(j-2, M)+1)) + (sigma(i, j) == sigma(i, mod(j, M)+1)));
    end
end
result = sum(H(:))/2;
for i = 1:N
    for j = 1:M
        result = result-H(i,j)*h;    
    end
end
end
% newtons method for PDE's
% written by fixedotter
n = 24; dx = 1/n;
I  = speye(n-1);
I_sq = speye((n-1)^2);
D2 = spdiags(kron(ones(n-1,1),[1 -2 1]),-1:1,n-1,n-1)/dx^2;
L  = kron(I,D2) + kron(D2, I);   % Laplacian = D_2^x + D_2^y

L_cub = kron(I_sq,D2) + kron(L,I); % lap with cube

[V,D] = eigs(L_cub,6,'sm'); % compute first 6 eigen pairs

diag(D)
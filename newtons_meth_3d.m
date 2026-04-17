% newtons method for PDE's
% written by fixedotter

% boring const
n = 10; dx = 1/n; ds = 0.1;

% basically just the grid
x = linspace(0,1,n+1)'; xs = x(2:end-1); % full and interior pts in interval 
[YS,XS,ZS] = meshgrid(xs,xs,xs); % grid of interior points

% defining the laplacian in 2d
I  = speye(n-1);
I_sq = speye((n-1)^2);
D2 = spdiags(kron(ones(n-1,1),[1 -2 1]),-1:1,n-1,n-1)/dx^2;

% Laplacian = D_2^x + D_2^y
L  = kron(I,D2) + kron(D2, I);

% laplacian w/ cube
L_cub = kron(I_sq,D2) + kron(L,I);

% compute first 6 eigen pairs
[V,D] = eigs(L_cub,6,'sm'); 
% diag(D)
% 6 smallest DIFFERENT eigenvalues are 3*pi^2, 6*pi^2, 9*pi^2, 11*pi^2, 12*pi^2, 14*pi^2
% since there's a linear combination, we basically get 3 of each of them after 3*pi^2, which is what we see in the diag(D) output

% function & it's jacobian with s to move along solution space?
F = @(s,u) L_cub*u + s*u + u.^3; % actual funct with s to move along the solutions
J = @(s,u) L_cub + spdiags(s + 3*u.^2, 0, (n-1)^3, (n-1)^3); % jacobian with s
f = @(x,y,z) sin(pi*x).*sin(pi*y).*sin(pi*z)*4; % init guess

% putting her into a column for easy matrix ops
u = reshape(f(XS,YS,ZS), (n-1)^3, 1);

% % do again with s & see if we can loop over the solutions to plot a grid
% for s = 0:ds:2*pi^2
%   chi = J(s,u)\F(s,u); % next newton step
%   u = u - chi; % update
%   if (norm(chi) < 1e-6) % stopping condition
%     break
%   end
% end

% just doing s = 1 for now
s = 1;
chi = J(s,u)\F(s,u); % next newton step
u = u - chi; % update
% if (norm(chi) < 1e-6) % stopping condition
%   disp("converged")
% else
%   disp("not converged")
% end

% bringing it back to 3d for plotting
u = reshape(u, n-1, n-1, n-1); % reshape function
Z = zeros(n+1,n+1,n+1); % boundary conditions
Z(2:end-1,2:end-1,2:end-1) = u; % feeding in the solution
[X,Y,Z] = meshgrid(x,x,x); % grid for plotting

% can we save this to a mf gif
gifFile = "neubergers_cube_w_newtonsMeth.gif";
for k = 0:n
  figure
  colormap(jet); % shading interp
  Z_slice = Z(:,:,k+1); % slicing the solution at z = k*ds
  zlim([-1,1]); % keeping the z axis fixed so it's less confusing
  surf(X(:,:,k+1),Y(:,:,k+1),Z_slice,'edgecolor','none'); % plotting the slice of the solution
  titleName = strcat('z = ', num2str(k*dx));
  title(titleName)
  view(42,-42);
  % exportgraphics(gcf, strcat('neubergers_cube_w_newtonsMeth_s=', num2str(k*ds), '.png')); % saving the figure
  exportgraphics(gcf,gifFile,Append=true); % saving to the gif
end


  % what chat said
  % figure('Visible','off','Name','gif');
  % colormap(jet); % shading interp
  % isosurface(X,Y,Z,k*ds); % plotting the isosurface of the solution
  % titleVar = strcat('s = ', num2str(k*ds));
  % title(titleVar)
  % view(42,-42);
  % saveas(gcf, strcat('wave_sol_s=', num2str(k*ds), '.png'));

% isosurface(X,Y,Z,0); % plotting the isosurface of the solution
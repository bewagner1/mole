clc
close all

addpath('../../src/matlab_octave')

% Parameters
k = 2;  % Order of accuracy
m = 40; % Number of nodes along x-axis
n = 40; % Number of nodes along y-axis
dx = 1 / (m - 1);
dy = 1 / (n - 1);
dc = [1; 1; 1; 1];
nc = [0; 0; 0; 0];

[X, Y] = genCurvGrid(n, m);
% [X, Y] = meshgrid(1:m, 1:n);

% Plot the physical grid
mesh(X, Y, zeros(n, m), 'Marker', '.', 'MarkerSize', 10)
view([0 90])
axis equal
set(gcf, 'Color', 'w')

% Staggered logical grid
[Xs, Ys] = meshgrid([1 1.5 : 1 : m-0.5 m], [1 1.5 : 1 : n-0.5 n]);
Cs = Xs.^2 + Ys.^2;
% Reshape the field so it can be multiplied by the operator later on
C_ = reshape(Cs', [], 1);

% Get 2D curvilinear mimetic gradient
G = grad2DCurv(k, Xs, Ys, m - 1, dx, n - 1, dy, dc, nc);

% Apply the operator to the field
TMP = G*C_;
Gx = TMP(1:(m+1)*(n));
Gy = TMP((m+1)*(n)+1:end);

% Reshape for visualization
Gx = reshape(Gx, m+1, n)';
Gy = reshape(Gy, m, n+1)';

CtoU = interpolCentersToFacesD1D(k, m - 1);
CtoU = kron(speye(n + 1), CtoU);
CtoV = interpolCentersToFacesD1D(k, n - 1);
CtoV = kron(CtoV, speye(m + 1));

Ux = CtoU * reshape(Xs', [], 1);
Uy = CtoU * reshape(Ys', [], 1);
Vx = CtoV * reshape(Xs', [], 1);
Vy = CtoV * reshape(Ys', [], 1);

Ux = reshape(Ux, m + 1, n)';
Uy = reshape(Uy, m + 1, n)';
Vx = reshape(Vx, m, n + 1)';
Vy = reshape(Vy, m, n + 1)';

% Plot results
figure
set(gcf, 'Color', 'w')
subplot(3, 1, 1)
surf(Xs, Ys, Cs, 'EdgeColor', 'none');
colorbar
xlabel('x')
ylabel('y')
title('C')
axis equal
view([0 90])
shading interp
subplot(3, 1, 2)
surf(Ux, Uy, Gx, 'EdgeColor', 'none');
colorbar
xlabel('x')
ylabel('y')
title('U')
axis equal
view([0 90])
shading interp
subplot(3, 1, 3)
surf(Vx, Vy, Gy, 'EdgeColor', 'none');
colorbar
xlabel('x')
ylabel('y')
title('V')
axis equal
view([0 90])
shading interp

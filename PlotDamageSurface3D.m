function []=A_PlotDamageSurface3D(alpha)

n_theta = 200; %nb of drawing points (different from the number of points
%used to compute the operators)
n_phi = n_theta;
[theta, phi] = meshgrid(linspace(0, pi, n_theta), linspace(0, 2*pi, n_phi));

n1 = sin(theta) .* cos(phi);
n2 = sin(theta) .* sin(phi);
n3 = cos(theta);

d_values = zeros(size(n1));

for i = 1:numel(n1)
    n = [n1(i); n2(i); n3(i)];
    N = n * n'; % tensor N_ij = n_i * n_j
    I = eye(3);

    % F_ij = N_ij - (1/3) * delta_ij
    F2 = N - (1/3) * I;

    % F_ijkl
    delta = @(i,j) double(i==j);
    F4 = zeros(3,3,3,3);
    for i1=1:3
        for j1=1:3
            for k1=1:3
                for l1=1:3
                    term1 = N(i1,j1) * N(k1,l1);
                    term2 = (1/7) * ( ...
                        delta(i1,j1)*N(k1,l1) + N(i1,j1)*delta(k1,l1) + ...
                        delta(i1,l1)*N(j1,k1) + N(i1,l1)*delta(j1,k1) + ...
                        delta(i1,k1)*N(j1,l1) + N(i1,k1)*delta(j1,l1));
                    term3 = (1/35) * ( ...
                        delta(i1,j1)*delta(k1,l1) + ...
                        delta(i1,l1)*delta(j1,k1) + ...
                        delta(i1,k1)*delta(j1,l1));
                    F4(i1,j1,k1,l1) = term1 - term2 + term3;
                end
            end
        end
    end

    % Construct Vd(n)
    Vd = zeros(21,1);
    Vd(1) = 1;
    Vh(2) = 1;

    % 2nd order terms
    Vd(3) = F2(1,1) - F2(3,3);
    Vd(4) = 2 * F2(1,2);
    Vd(5) = 2 * F2(1,3);
    Vd(6) = F2(2,2) - F2(3,3);
    Vd(7) = 2 * F2(2,3);

    Vh(8:12) = Vd(3:7); % Same components than Vd

    % 4th order terms (indices for F_ijkl)
    getF = @(i,j,k,l) F4(i,j,k,l);

    Vd(13) = getF(1,1,1,1) + getF(3,3,3,3) - 2*getF(1,1,3,3) - 4*getF(1,3,1,3);
    Vd(14) = 2*getF(3,3,3,3) + 2*getF(1,1,2,2) - 2*getF(1,1,3,3) - 2*getF(2,2,3,3) ...
        + 4*getF(1,2,1,2) - 4*getF(1,3,1,3) - 4*getF(2,3,2,3);
    Vd(15) = getF(2,2,2,2) + getF(3,3,3,3) - 2*getF(2,2,3,3) - 4*getF(2,3,2,3);
    Vd(16) = 4*(getF(1,1,1,2) - getF(3,3,1,2) - 2*getF(1,2,2,3) - 2*getF(1,3,2,3));
    Vd(17) = 4*(getF(1,1,1,3) - getF(3,3,1,3));
    Vd(18) = 4*(getF(1,1,2,3) - getF(3,3,2,3) + 2*getF(1,2,1,3));
    Vd(19) = 4*(getF(2,2,1,2) - getF(3,3,1,2) - 2*getF(1,2,2,3) - 2*getF(1,3,2,3));
    Vd(20) = 4*(getF(2,2,1,3) - getF(3,3,1,3));
    Vd(21) = 4*(getF(2,2,2,3) - getF(3,3,2,3));

    d_values(i) = dot(Vd, alpha);
end

R_d = abs(d_values);

% Cartesian coordinates of the surface
X_d = R_d .* n1;
Y_d = R_d .* n2;
Z_d = R_d .* n3;

figure

% d(n) surface
surf(X_d, Y_d, Z_d, d_values, 'EdgeColor', 'none');
colormap(gca, jet);
title('Damage orientation function d(\theta,\phi)');
axis equal;
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);  % Boîte fixe
colorbar; view(3);
xlabel('x'); ylabel('y'); zlabel('z');
box
lighting gouraud;
camlight headlight;
material shiny;

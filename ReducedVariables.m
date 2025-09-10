function []=ReducedVariables(SAVE_Alphas)

%PCA with automatic selection based on error criterion
Error=0.01;%tolerance
P = 100-Error; 
X = SAVE_Alphas;  %Original data (N x D), not centred

%Covariance
C = (X' * X) / size(X,1);  % (D x D)

% SVD
[U, S, ~] = svd(C);

% Compute the energy
eigenvalues = diag(S);  
explained = 100 * eigenvalues / sum(eigenvalues); 
cumulative = cumsum(explained); 

% Find the minimum number of components M such that variance ≥ P%.
M = find(cumulative >= P, 1, 'first');

%Projection on the first M principal components
W = U(:, 1:M);        % (D x M)
X_reduced = X * W;    % (N x M)

%Plot reduced internal variables
figure;
hold on;
colors = lines(M);
for i = 1:M
    plot(1:size(X,1), -X_reduced(:,i), 'LineWidth', 2, 'Color', colors(i,:));
end
%Note: beta variables are plotted arbitrarily with a minus sign to match the
%evolution sign of alphas, even though the sign of beta does not
%have any influence when used in a reduced database

xlabel('n','Fontsize',16);
ylabel('\beta_i','Fontsize',16);
title('Reduced internal variable(s) \beta_i by PCA');
grid on;
box



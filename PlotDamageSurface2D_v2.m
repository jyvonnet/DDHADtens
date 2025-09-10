function []=A_PlotDamageSurface2D_v2(alpha)




theta = linspace(0, 2*pi, 360);   
d_theta = zeros(size(theta));     

for i = 1:length(theta)
    th = theta(i);
    n1 = cos(th);
    n2 = sin(th);

    F11 = n1^2 - 0.5;
    F22 = n2^2 - 0.5;
    F12 = n1 * n2;

    F1111 = n1^4 - n1^2 +  1/8;
    F1122 = n1^2 * n2^2 - 1/8;
    F2222 = n2^4 - n2^2 + 1/8;
    F1212 = n1^2 * n2^2 - 1/8;

    F1112 = n1^3 * n2 - (1/2)*n1*n2;
    F2212 = n1*n2^3  - (1/2)*n1*n2;

%    Construction de Vd(θ) selon Eq. (21)
    Vd = [
        1;
        0;
        F11 - F22;
        2 * F12;
        F1111 - 2*F1122 + F2222 - 4*F1212;
        4 * (F1112 - F2212)
    ];

    d_theta(i) = Vd' * alpha';
end

figure
polarplot(theta, d_theta, 'r-', 'LineWidth', 2);
title('Damage orientation function d(\theta)');
rlim([0 1])
ax = gca;
ax.ThetaZeroLocation = 'right';
ax.ThetaDir = 'counterclockwise';
grid on;



%DDHADtens 2D v1
%J. Yvonnet July 1st, 2025
%Gustave Eiffel University
%Input : a collection of 3x3 elastic tensors
%Output: associated collection of damage variables (harmonic analysis)
%6x1 vector
%plots of damage variable evolution and damage surfaces (2D)
%-----------------------------------------------------------
%The database is associated with the numerical examples in 
%[1] J. Yvonnet, Q.-C. He, P. Li, A data-driven harmonic approach to
% constructing anisotropic damage models with a minimum number of internal
% variables, Journal of the Mechanics and Physics of Solids,
% 162:104828, 2022. 
%https://doi.org/10.1016/j.jmps.2022.104828

clear all
close all

%1) Load database
%=================
[CC_save,CC_0]=LoadDatabase_DDHAD2D;

%2) Extract internal variables
%==================================
Npoints=50; %number of points for discretizing the unit circle
Dthet=2*pi/(Npoints-1); %angular increment
thet =-pi:Dthet:pi;

[eta0,kappa0] = CaculEtaKappa0_2D(CC_save,thet); 
[C_tilde] = PrecalculOPE_2D(eta0,kappa0,thet);
[SAVE_Alphas] = CalculVarDirect_2D(C_tilde,CC_0,CC_save);

hold on
plot(SAVE_Alphas,'LineWidth',2) %show internal variables
xlabel('n','Fontsize',16)
ylabel('\alpha_i','Fontsize',16)
title('Internal variables \alpha_i');
box
grid on


%3) Compare initial tensor components and reconstruction 
%from internal variables
%=======================================================
C_vec0=[CC_0(1,1) CC_0(1,2) CC_0(1,3) CC_0(2,2) CC_0(2,3) CC_0(3,3)]';

for i=1:size(CC_save,3)
    CC=CC_save(:,:,i);
    alpha=SAVE_Alphas(i,:);
    C_vec(:,i)=[CC(1,1) CC(1,2) CC(1,3) CC(2,2) CC(2,3) CC(3,3)]';%initial tensor
    CC_save_REC(:,i)=C_vec0-C_tilde*alpha';%reconstructed tensor from internal variables
end
figure
hold on
h = plot(C_vec', 'LineWidth', 2);
step = 5;
colors = vertcat(h.Color); 
idx = 1:step:size(CC_save_REC, 2);
set(plot(idx,CC_save_REC(:,idx)', 'o', 'LineStyle', 'none', ...
    'LineWidth', 2), ...
    {'MarkerEdgeColor'}, num2cell(colors, 2));

xlabel('n','Fontsize',16)
ylabel('C_{ijkl}','Fontsize',16)
title('Elastic tensor components');
box on
grid on
dummy_line = plot(NaN, NaN, '-', 'Color', 'k', 'LineWidth', 2); 
dummy_marker = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'k', 'LineStyle', 'none'); 
legend([dummy_line, dummy_marker], {'Data', 'Reconstructed'}, 'Location', 'best');

%4) Show damage surface
%===========================
PlotDamageSurface2D_v2(alpha)

%5) Reduced internal variables with PCA
%=====================================
%note : this reduction is here performed only along one
%loading. If more samples of elastic tensors arising from
%a larger collection of loads is available, the reduction can be
%performed on the whole set of data (collection of all damage
%variables alpha history)
ReducedVariables(SAVE_Alphas);








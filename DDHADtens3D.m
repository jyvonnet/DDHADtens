clear all
close all
%DDHADtens 3D v1
%J. Yvonnet July 1st, 2025
%Gustave Eiffel University
%Input : a collection of 6x6 elastic tensors
%Output: associated collection of damage variables (harmonic analysis)
%21x1 vector
%plots of damage variable evolution and damage surfaces (2D)
%-----------------------------------------------------------
%The database is associated with the numerical examples in 
%[2] J. Yvonnet, Qi-Chang He, Microstructure-based machine learning of damage
% models including anisotropy, irreversibility and evolution, 
% Journal of the Mechanics and Physics of Solids, 106160, 2025.

%1) Load database
%=================
[CC_save,CC_0]=LoadDatabase_DDHAD3D;


%2) Extract internal variables
%==================================
NbPts=50; %number of points for discretizing the unit sphere
Delll=pi/(NbPts); %angular increment
thet=0:Delll:pi; phi=0:Delll:2*pi;

[eta0,kappa0] = CaculEtaKappa0_3D(CC_0,thet,phi);
[C_tilde] = PrecalculOPE_3D(eta0,kappa0,thet,phi);
[SAVE_Alphas] = CalculVarDirect_3D(C_tilde,CC_0,CC_save);
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
C_vec0=[CC_0(1,1) CC_0(1,2) CC_0(1,3) CC_0(1,4) CC_0(1,5) CC_0(1,6) ...
    CC_0(2,2) CC_0(2,3) CC_0(2,4) CC_0(2,5) CC_0(2,6)  ...
    CC_0(3,3) CC_0(3,4) CC_0(3,5) CC_0(3,6) CC_0(4,4) CC_0(4,5) CC_0(4,6) ...
    CC_0(5,5) CC_0(5,6) CC_0(6,6)]';

for i=1:size(CC_save,3)
    CC=CC_save(:,:,i);
    alpha=SAVE_Alphas(i,:);
    C_vec(:,i)=[CC(1,1) CC(1,2) CC(1,3) CC(1,4) CC(1,5) CC(1,6) ...
    CC(2,2) CC(2,3) CC(2,4) CC(2,5) CC(2,6)...
    CC(3,3) CC(3,4) CC(3,5) CC(3,6) ...
    CC(4,4) CC(4,5) CC(4,6) ...
    CC(5,5) CC(5,6) CC(6,6)]';%initial tensor
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

%4)Show evolution of damage
%==========================
alpha=SAVE_Alphas(size(SAVE_Alphas,1),:); %last damage state
PlotDamageSurface3D(alpha);

%5) Extract reduced variables with PCA
%=====================================
%note : this reduction is here performed only along one
%loading. If more samples of elastic tensors arising from
%a larger collection of loads is available, the reduction can be
%performed on the whole set of data (collection of all damage
%variables alpha history)
ReducedVariables(SAVE_Alphas);








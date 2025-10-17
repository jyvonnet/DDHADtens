function [CC_save,CC_0]=A_LoadDatabase_DDHAD3D()

fprintf('==============================\n');
fprintf('1: RVE A: RVE with brittle matrix, hard inclusion, traction\n');
fprintf('2: RVE B: RVE with weak layer perpendicular to x-axis, traction\n');
fprintf('3: RVE C: Porous RVE, compression\n');
fprintf('4: RVE D: Fiber in brittle matrix directed along x, traction\n');
fprintf('==============================\n');
Type=input('Choose database (1-4): ');

if Type==1
    load('RVEA_3D');
elseif Type==2
    load('RVEB_3D');
elseif Type==3
    load('RVEC_3D');
elseif Type==4
    load('RVED_3D');
end
CC_0=CC_save2(:,:,1); %undamaged elastic tensor (6x6 matrix)
CC_save=CC_save2; %history of elastic tensors (6x6 matrices)







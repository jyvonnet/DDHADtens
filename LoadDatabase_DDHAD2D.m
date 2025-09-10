function [CC_save,CC_0]=A_LoadDatabase_DDHAD2D()

'=============================='
'1: RVE A: RVE with vertical weak layer, loading along x'
'2: RVE B: RVE with circular hole, loading along x'
'3: RVE C: RVE with circular inclusion, loading along x'
'4: RVE D: concrete microstructure (image), loading along x'
'5: RVE E: RVE with 45° weak layer, loading along x'
'=============================='
Type=input('Choose database')


if Type==1
    load('RVEA_2D');
elseif Type==2
    load('RVEB_2D')
elseif Type==3
    load('RVEC_2D')
elseif Type==4
    load('RVED_2D')
elseif Type==5
    load('RVEE_2D')
end


CC_0=CC_save2(:,:,1); %undamaged elastic tensor (3x3 matrix)
CC_save=CC_save2; %history of elastic tensors (3x3 matrices)







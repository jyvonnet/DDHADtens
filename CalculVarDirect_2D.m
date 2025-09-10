function [SAVE_Alphas] = A_CalculVarDirect_2D(C_tilde,CC_0,CC_save)

%given C_tilde,CC_0, compute \alpha for each C in CC_save

C_0_vec=[CC_0(1,1) CC_0(1,2) CC_0(1,3) CC_0(2,2) CC_0(2,3) CC_0(3,3)]';

K=C_tilde'*C_tilde;

SAVE_Alphas=zeros(size(CC_save,3),6);
for i=1:size(CC_save,3)
    Cx=CC_save(:,:,i);
    C_x_vec=[Cx(1,1) Cx(1,2) Cx(1,3) Cx(2,2) Cx(2,3) Cx(3,3) ]';
    G=C_tilde'*(C_0_vec-C_x_vec);
    alpha_i=K\G;
    SAVE_Alphas(i,:)=alpha_i';
end




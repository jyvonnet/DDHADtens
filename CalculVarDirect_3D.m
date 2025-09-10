function [SAVE_Alphas] = A_CalculVarDirect_3D(C_tilde,CC_0,CC_save)

%given C_tilde,CC_0, compute \alpha for each C in CC_save

C_0_vec=[CC_0(1,1) CC_0(1,2) CC_0(1,3) CC_0(1,4) CC_0(1,5) CC_0(1,6) CC_0(2,2) CC_0(2,3) CC_0(2,4) CC_0(2,5) CC_0(2,6) ...
    CC_0(3,3) CC_0(3,4) CC_0(3,5) CC_0(3,6) CC_0(4,4) CC_0(4,5) CC_0(4,6) CC_0(5,5) CC_0(5,6) CC_0(6,6)]';

K=C_tilde'*C_tilde;

SAVE_Alphas=zeros(size(CC_save,3),21);
for i=1:size(CC_save,3)
    Cx=CC_save(:,:,i);
    C_x_vec=[Cx(1,1) Cx(1,2) Cx(1,3) Cx(1,4) Cx(1,5) Cx(1,6) Cx(2,2) Cx(2,3) Cx(2,4) Cx(2,5) Cx(2,6) ...
        Cx(3,3) Cx(3,4) Cx(3,5) Cx(3,6) Cx(4,4) Cx(4,5) Cx(4,6) Cx(5,5) Cx(5,6) Cx(6,6)]';
    G=C_tilde'*(C_0_vec-C_x_vec);
    alpha_i=K\G;
    SAVE_Alphas(i,:)=alpha_i';
end




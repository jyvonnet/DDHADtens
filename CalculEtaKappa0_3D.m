function [eta0,kappa0] = A_CalculEtaKappa0_3D(MatC,thet,phi)

%compute undamaged direction-dependent elastic parameters
%eta^0 = n_i n_j n_k n_l C_ijkl^0
%kappa^0 = \delta_ij n_k n_l C_ijkl^0

C1111=MatC(1,1);
C1122=MatC(1,2);
C1133=MatC(1,3);
C1112=MatC(1,4);
C1113=MatC(1,5);
C1123=MatC(1,6);

C2222=MatC(2,2);
C2233=MatC(2,3);
C2212=MatC(2,4);
C2213=MatC(2,5);
C2223=MatC(2,6);

C3333=MatC(3,3);
C3312=MatC(3,4);
C3313=MatC(3,5);
C3323=MatC(3,6);

C1212=MatC(4,4);
C1213=MatC(4,5);
C1223=MatC(4,6);

C1313=MatC(5,5);
C1323=MatC(5,6);

C2323=MatC(6,6);

NT=length(thet)-1;
MP=length(phi)-1;
MatTeta=zeros(NT,MP);
MatPhi=zeros(NT,MP);
eta=zeros(NT,MP);
kappa=zeros(NT,MP);

for i=1:length(thet)-1
    
    for j=1:length(phi)-1
        TT=thet(i);
        PP=phi(j);
        
        MatTeta(i,j)=TT;
        MatPhi(i,j)=PP;
        
        n1=sin(TT)*cos(PP);
        n2=sin(TT)*sin(PP);
        n3=cos(TT);

        eta(i,j) = C1111*n1^4 +C2222*n2^4+C3333*n3^4+...
            2*C1122*n1^2*n2^2+2*C1133*n1^2*n3^2+2*C2233*n2^2*n3^2+...
            4*C1212*n1^2*n2^2+4*C1313*n1^2*n3^2+4*C2323*n2^2*n3^2+...
            4*C1112*n1^3*n2+4*C1113*n1^3*n3+4*C2212*n2^3*n1+4*C3313*n1*n3^3+...
            4*C1123*n1^2*n2*n3+4*C2213*n2^2*n1*n3+4*C3312*n1*n2*n3^2+4*C2223*n2^3*n3+4*C3323*n2*n3^3+...
            8*C1213*n1^2*n2*n3+  8*C1223*n1*n2^2*n3+8*C1323*n1*n2*n3^2;
                       
        kappa(i,j)=C1111*n1^2+C2222*n2^2+C3333*n3^2+...
            C1122*n1^2+ C1122*n2^2+...
            C1133*n3^2+ C1133*n1^2 +...
            +C2233*n3^2++C2233*n2^2+...
            +2*C1112*n1*n2+2*C1113*n1*n3+2*C1123*n2*n3+...
            +2*C2212*n1*n2+2*C2213*n1*n3+2*C2223*n2*n3+...
            +2*C3312*n1*n2+2*C3313*n1*n3+2*C3323*n2*n3;
   
    end
end
eta0=eta;
kappa0=kappa;


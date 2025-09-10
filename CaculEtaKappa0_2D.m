function [eta0,kappa0] = A_CaculEtaKappa0_2D(CC_save,thet)

C1111=CC_save(1,1);
C2222=CC_save(2,2);
C1212=CC_save(3,3);
C2212=CC_save(2,3);
C1112=CC_save(1,3);
C1122=CC_save(1,2);

NT=length(thet)-1;
MatTeta=zeros(NT,1);

eta=zeros(NT,1);
kappa=zeros(NT,1);

for i=1:length(thet)-1

    TT=thet(i);

    MatTeta(i,1)=TT;

    n1=cos(TT);
    n2=sin(TT);

    eta(i,1) = C1111*n1^4 +C2222*n2^4+...
        2*C1122*n1^2*n2^2+...
        4*C1212*n1^2*n2^2+...
        4*C1112*n1^3*n2+4*C2212*n2^3*n1;

    kappa(i,1)=C1111*n1^2+C2222*n2^2+...
        C1122*n1^2+ C1122*n2^2+...
        +2*C1112*n1*n2+...
        +2*C2212*n1*n2;

end

eta0=eta;
kappa0=kappa;



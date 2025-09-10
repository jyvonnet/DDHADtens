function [C_tilde] = A_PrecalculOPE_2D(eta0,kappa0,thet)


Dthet=2*pi/(length(thet)-1);

Somme_vecud=zeros(6,1);
Somme_vecvd=zeros(6,1);
Somme_vecV11d=zeros(6,1);
Somme_vecV12d=zeros(6,1);
Somme_vecZ1111d=zeros(6,1);
Somme_vecZ1112d=zeros(6,1);

for i=1:length(thet)-1
    
    TT=thet(i);
   
    n1=cos(TT);
    n2=sin(TT);
   
    F11=(n1^2-1/2);
    F22=(n2^2-1/2);
    F12=(n1*n2);
    
    F1111 = (   n1^4-n1^2+1/8   );
    F1122 = n1^2*n2^2 - (1/8);
    F1112 = (    n1^3*n2 - (1/2)*n1*n2    );
    
    F2222 = n2^4-n2^2+1/8;
    F2212 = n1*n2^3 -  (1/2)*(n1*n2);
    
    F1212 = n1^2*n2^2 - (1/8);
    
    VecD=[1
        0
        F11-F22
        2*F12
        F1111-2*F1122+F2222-4*F1212
        4*F1112-4*F2212];
    
    VecH=[0
        1
        (F11-F22);
        2*F12;
        0
        0];
    
    Somme_vecud=Somme_vecud+eta0(i)*VecD*Dthet;
    Somme_vecvd=Somme_vecvd+kappa0(i)*VecH*Dthet;
    Somme_vecV11d=Somme_vecV11d+F11*eta0(i)*VecD*Dthet;
    Somme_vecV12d=Somme_vecV12d+F12*eta0(i)*VecD*Dthet;
    Somme_vecZ1111d=Somme_vecZ1111d+F1111*eta0(i)*VecD*Dthet;
    Somme_vecZ1112d=Somme_vecZ1112d+F1112*eta0(i)*VecD*Dthet;
 
end

vecud =(1/(2*pi))*Somme_vecud;
vecvd = (1/(2*pi))*Somme_vecvd;
vecV11d = (2/(pi))*Somme_vecV11d;
vecV12d = (2/(pi))*Somme_vecV12d;
vecZ1111d = (8/(pi))*Somme_vecZ1111d;
vecZ1112d = (8/(pi))*Somme_vecZ1112d;

C1111tt=vecud+vecV11d+vecZ1111d;
C1122tt=vecvd-vecud-vecZ1111d;
C1112tt=0.5*vecV12d+vecZ1112d;
C2222tt=vecud-vecV11d+vecZ1111d;
C2212tt=0.5*vecV12d-vecZ1112d;
C1212tt=vecud-0.5*vecvd-vecZ1111d;

C_tilde=[C1111tt C1122tt C1112tt C2222tt C2212tt C1212tt]';




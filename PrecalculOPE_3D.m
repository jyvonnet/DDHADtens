function [C_tilde] = A_PrecalculOPE_3D(eta0,kappa0,thet,phi)

%compute C_tilde operator such that
%C(t) = C^0 - \sum_k \alpha_k C_tilde^k

Dthet=pi/(length(thet)-1);
Dphi=2*pi/(length(phi)-1);

Somme_vecad=zeros(21,1);
Somme_vecbh=zeros(21,1);

Somme_vecU11d=zeros(21,1);
Somme_vecU12d=zeros(21,1);
Somme_vecU13d=zeros(21,1);
Somme_vecU22d=zeros(21,1);
Somme_vecU23d=zeros(21,1);

Somme_vecV11h=zeros(21,1);
Somme_vecV12h=zeros(21,1);
Somme_vecV13h=zeros(21,1);
Somme_vecV22h=zeros(21,1);
Somme_vecV23h=zeros(21,1);

Somme_vecZ1111d=zeros(21,1);
Somme_vecZ1122d=zeros(21,1);
Somme_vecZ2222d=zeros(21,1);
Somme_vecZ1112d=zeros(21,1);
Somme_vecZ1113d=zeros(21,1);
Somme_vecZ1123d=zeros(21,1);
Somme_vecZ2212d=zeros(21,1);
Somme_vecZ2213d=zeros(21,1);
Somme_vecZ2223d=zeros(21,1);


for i=1:length(thet)-1

    for j=1:length(phi)-1
        TT=thet(i);
        PP=phi(j);


        n1=sin(TT)*cos(PP);
        n2=sin(TT)*sin(PP);
        n3=cos(TT);

        F11=(n1^2-1/3);
        F12=(n1*n2);
        F13=(n1*n3);
        F22=(n2^2-1/3);
        F23=(n2*n3);
        F33=(n3^2-1/3);


        F1111 = n1^4-(6/7)*n1^2+3/35;
        F1122 = n1^2*n2^2 - (1/7)*(n1^2+n2^2)+1/35;
        F1133 = n1^2*n3^2 - (1/7)*(n1^2+n3^2)+1/35;
        F1112 = n1^3*n2 - (3/7)*n1*n2;
        F1113 = n1^3*n3 - (3/7)*n1*n3;
        F1123 = n1^2*n2*n3 - (1/7)*n2*n3;

        F2222 = n2^4-(6/7)*n2^2+3/35;
        F2233 = n2^2*n3^2 - (1/7)*(n2^2+n3^2)+1/35;
        F2212 = n1*n2^3 -  (3/7)*(n1*n2);
        F2213 = n1*n2^2*n3 - (1/7)*(n1*n3); %%
        F2223 = n2^3*n3 - (3/7)*n2*n3;

        F3333 = n3^4-(6/7)*n3^2+3/35;
        F3312 = n1*n2*n3^2 - (1/7)*n1*n2;
        F3313 = n1*n3^3 - (3/7)*n1*n3; %%
        F3323 = n2*n3^3 - (3/7)*n2*n3;

        F1212 = n1^2*n2^2 - (1/7)*(n1^2+n2^2)+1/35;
        F1213 = n1^2*n2*n3 - (1/7)*n2*n3; %%
        F1223 = n1*n2^2*n3 - (1/7)*n1*n3; %%

        F1313 = n1^2*n3^2 - (1/7)*(n1^2+n3^2)+1/35;
        F1323 = n1*n2*n3^2 - (1/7)*n1*n2;

        F2323 = n2^2*n3^2 - (1/7)*(n2^2+n3^2)+1/35;

        VecD=[1
            0
            F11-F33
            2*F12
            2*F13
            F22-F33
            2*F23
            0
            0
            0
            0
            0
            F1111+F3333-2*F1133-4*F1313
            2*F3333+2*F1122-2*F1133-2*F2233+4*F1212-4*F1313-4*F2323
            F2222+F3333-2*F2233-4*F2323
            4*F1112-4*F3312-8*F1223-8*F1323
            4*F1113-4*F3313
            4*F1123-4*F3323+8*F1213
            4*F2212-4*F3312-8*F1223-8*F1323
            -4*F3313+4*F2213
            4*F2223-4*F3323];

        VecH=[0
            1
            0
            0
            0
            0
            0
            F11-F33
            2*F12
            2*F13
            F22-F33
            2*F23
            0
            0
            0
            0
            0
            0
            0
            0
            0];

        Somme_vecad=Somme_vecad+eta0(i,j)*VecD*Dthet*Dphi*sin(TT);
        Somme_vecbh=Somme_vecbh+kappa0(i,j)*VecH*Dthet*Dphi*sin(TT);

        Somme_vecU11d=Somme_vecU11d+eta0(i,j)*VecD*F11*Dthet*Dphi*sin(TT);
        Somme_vecU12d=Somme_vecU12d+eta0(i,j)*VecD*F12*Dthet*Dphi*sin(TT);
        Somme_vecU13d=Somme_vecU13d+eta0(i,j)*VecD*F13*Dthet*Dphi*sin(TT);
        Somme_vecU22d=Somme_vecU22d+eta0(i,j)*VecD*F22*Dthet*Dphi*sin(TT);
        Somme_vecU23d=Somme_vecU23d+eta0(i,j)*VecD*F23*Dthet*Dphi*sin(TT);

        Somme_vecV11h=Somme_vecV11h+kappa0(i,j)*VecH*F11*Dthet*Dphi*sin(TT);
        Somme_vecV12h=Somme_vecV12h+kappa0(i,j)*VecH*F12*Dthet*Dphi*sin(TT);
        Somme_vecV13h=Somme_vecV13h+kappa0(i,j)*VecH*F13*Dthet*Dphi*sin(TT);
        Somme_vecV22h=Somme_vecV22h+kappa0(i,j)*VecH*F22*Dthet*Dphi*sin(TT);
        Somme_vecV23h=Somme_vecV23h+kappa0(i,j)*VecH*F23*Dthet*Dphi*sin(TT);


        Somme_vecZ1111d=Somme_vecZ1111d+eta0(i,j)*VecD*F1111*Dthet*Dphi*sin(TT);
        Somme_vecZ1122d=Somme_vecZ1122d+eta0(i,j)*VecD*F1122*Dthet*Dphi*sin(TT);
        Somme_vecZ2222d=Somme_vecZ2222d+eta0(i,j)*VecD*F2222*Dthet*Dphi*sin(TT);
        Somme_vecZ1112d=Somme_vecZ1112d+eta0(i,j)*VecD*F1112*Dthet*Dphi*sin(TT);
        Somme_vecZ1113d=Somme_vecZ1113d+eta0(i,j)*VecD*F1113*Dthet*Dphi*sin(TT);
        Somme_vecZ1123d=Somme_vecZ1123d+eta0(i,j)*VecD*F1123*Dthet*Dphi*sin(TT);
        Somme_vecZ2212d=Somme_vecZ2212d+eta0(i,j)*VecD*F2212*Dthet*Dphi*sin(TT);
        Somme_vecZ2213d=Somme_vecZ2213d+eta0(i,j)*VecD*F2213*Dthet*Dphi*sin(TT);
        Somme_vecZ2223d=Somme_vecZ2223d+eta0(i,j)*VecD*F2223*Dthet*Dphi*sin(TT);





    end
end


Somme_vecad=(1/(4*pi))*Somme_vecad;
Somme_vecbh=(1/(4*pi))*Somme_vecbh;

Somme_vecU11d=(15/(8*pi))*Somme_vecU11d;
Somme_vecU12d=(15/(8*pi))*Somme_vecU12d;
Somme_vecU13d=(15/(8*pi))*Somme_vecU13d;
Somme_vecU22d=(15/(8*pi))*Somme_vecU22d;
Somme_vecU23d=(15/(8*pi))*Somme_vecU23d;

Somme_vecV11h=(15/(8*pi))*Somme_vecV11h;
Somme_vecV12h=(15/(8*pi))*Somme_vecV12h;
Somme_vecV13h=(15/(8*pi))*Somme_vecV13h;
Somme_vecV22h=(15/(8*pi))*Somme_vecV22h;
Somme_vecV23h=(15/(8*pi))*Somme_vecV23h;


Somme_vecZ1111d=(315/(32*pi))*Somme_vecZ1111d;
Somme_vecZ1122d=(315/(32*pi))*Somme_vecZ1122d;
Somme_vecZ2222d=(315/(32*pi))*Somme_vecZ2222d;
Somme_vecZ1112d=(315/(32*pi))*Somme_vecZ1112d;
Somme_vecZ1113d=(315/(32*pi))*Somme_vecZ1113d;
Somme_vecZ1123d=(315/(32*pi))*Somme_vecZ1123d;
Somme_vecZ2212d=(315/(32*pi))*Somme_vecZ2212d;
Somme_vecZ2213d=(315/(32*pi))*Somme_vecZ2213d;
Somme_vecZ2223d=(315/(32*pi))*Somme_vecZ2223d;




alpha=(1/2)*(Somme_vecbh-Somme_vecad);
beta=(1/4)*(3*Somme_vecad-Somme_vecbh);

AA11d=Somme_vecV11h-Somme_vecU11d;
AA22d=Somme_vecV22h-Somme_vecU22d;
AA12d=Somme_vecV12h-Somme_vecU12d;
AA13d=Somme_vecV13h-Somme_vecU13d;
AA23d=Somme_vecV23h-Somme_vecU23d;

BB11d=(1/4)*(3*Somme_vecU11d-2*Somme_vecV11h);
BB22d=(1/4)*(3*Somme_vecU22d-2*Somme_vecV22h);
BB12d=(1/4)*(3*Somme_vecU12d-2*Somme_vecV12h);
BB13d=(1/4)*(3*Somme_vecU13d-2*Somme_vecV13h);
BB23d=(1/4)*(3*Somme_vecU23d-2*Somme_vecV23h);



C1111tt=alpha + 2*beta+2*AA11d+4*BB11d+Somme_vecZ1111d;
C1122tt=alpha +(AA11d+AA22d)+Somme_vecZ1122d;
C1133tt=alpha -AA22d-Somme_vecZ1111d-Somme_vecZ1122d;
C1112tt= AA12d+2*BB12d+Somme_vecZ1112d;
C1113tt= AA13d+2*BB13d+Somme_vecZ1113d;
C1123tt= AA23d+Somme_vecZ1123d;


C2222tt=alpha + 2*beta+2*AA22d+4*BB22d+Somme_vecZ2222d;
C2233tt=alpha -AA11d-Somme_vecZ1122d-Somme_vecZ2222d;
C2212tt= AA12d+2*BB12d+Somme_vecZ2212d;
C2213tt= AA13d+Somme_vecZ2213d;
C2223tt= AA23d+2*BB23d+Somme_vecZ2223d;

C3333tt=alpha + 2*beta-2*(AA11d+AA22d)-4*(BB11d+BB22d)+Somme_vecZ1111d+Somme_vecZ2222d+2*Somme_vecZ1122d;
C3312tt= AA12d-Somme_vecZ1112d-Somme_vecZ2212d;
C3313tt= AA13d+2*BB13d-Somme_vecZ1113d-Somme_vecZ2213d;
C3323tt= AA23d+2*BB23d-Somme_vecZ1123d-Somme_vecZ2223d;

C1212tt=beta+(BB11d+BB22d)+Somme_vecZ1122d;
C1213tt= BB23d+Somme_vecZ1123d;
C1223tt=BB13d+Somme_vecZ2213d;

C1313tt=beta-BB22d-Somme_vecZ1111d-Somme_vecZ1122d;
C1323tt=BB12d-Somme_vecZ1112d-Somme_vecZ2212d;

C2323tt=beta-BB11d-Somme_vecZ1122d-Somme_vecZ2222d;

C_tilde=[C1111tt C1122tt C1133tt C1112tt C1113tt C1123tt C2222tt C2233tt C2212tt C2213tt C2223tt C3333tt C3312tt C3313tt C3323tt C1212tt C1213tt C1223tt C1313tt C1323tt C2323tt]';



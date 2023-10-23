% Méthode de Trajectoire Centrale via une nouvelle fonction noyau pour l'optimisation convexe (Linéaire) 
%clc
%clear all
tic
%Les données
%Le paramètre de mise à jour
theta=0.99;
%le paramètre beta pour le pas déplacement
beta=0.99;
%Le paramètre de précision
 epsilon=10^(-4);
% Exemple de Netlib
load("lp_scagr7.mat")
A = Problem.A; 
a=full(A);
b = Problem.b;
b=full(b);
c = Problem.aux.c;
c=full(c);
[m,n]=size(A);
%Le point initial est trouvé par l'exécution de fichier PrimalDualPL.m
% x=zeros(n,1);
% z=zeros(n,1);
 % y=zeros(m,1);
%le paramètre de la fonction noyau
 p=log(n)/2;
 %Le paramètre de limite
tau=sqrt(n);
 %le paramètre barrière
 mu=x'*z/n;
%le nombre des itérations externes dans l'algorithme
kext=0;
%le nombre des itération internes dans l'algorithme
kint=0;
%La fonction noyau psi=inline('(t^2-1)/2+log(2/(1+t))-(exp(1/(1+t)-0.5)/(2*p))-(1/(2*p))');
 psi=inline('(t^2-1)/2+log(2/(1+t))-(exp(1/(1+t)-0.5)/(2*p))-(1/(2*p))');
psiv=zeros(n,1);
%La dérivée première de la fonction noyau psi 
psi1=inline('t-(1/(1+t))-(2*(exp(1/(1+t)-0.5)/(2*p))/((1+t)^2))');
psi1v=zeros(n,1);
 while x'*z>epsilon
     %itération externe
     mu=(1-theta)*mu;
     %le vecteur réduit
     v=sqrt(x.*z/mu);
     %La fonction de proximité phi(v)
     for i=1:n
         psiv(i)=psi(p,v(i));
     end
     phi=sum(psiv);
   
     while phi>tau
       %itération externe
       %%la résolution du système de Newton
       %le gradient de phi(v)
       for i=1:n
           %psi1v(i)=psi1(p,q,v(i));
           psi1v(i)=psi1(p,v(i));
       end
       r=-psi1v;
       s=zeros(m,1);
       t=zeros(n,1);
       %la décomposition de Cholesky pour trouver deltay%
       a1=1/mu*a*(diag(v)^(-1))*diag(x);
       M=a1*a1';
       R=chol(M);
       s1=s-a1*(r+t);
       u=R'\s1;
       deltay=R\u;
       %
       dx=a1'*deltay+(r+t);
       dz=-a1'*deltay-t;
       %la direction de Newton classique
       deltax=x.*dx./v;
       deltaz=z.*dz./v;
       %le pas de déplacement
       delta=1/2*norm(psi1v,2);
  %le pas de déplacement 
      alphax=depl(x,deltax);
      alphaz=depl(z,deltaz);
     alpha=beta*min(alphax,alphaz);
       %le nouvel itéré
       x=x+alpha*deltax;
       y=y+alpha*deltay;
       z=z+alpha*deltaz;
       %le nombre d'itération interne
       kint=kint+1;
       %le vecteur réduit
       v=sqrt(x.*z/mu);
       %La fonction de proximité phi(v)
       for i=1:n
         psiv(i)=psi(p,v(i));
       end
       phi=sum(psiv);
     end
     kext=kext+1;
 end
 toc
 primal=c'*x;
 dual=b'*y;
 saut=x'*z;
% fprintf('La valeur de saut est', saut);
  %fprintf('La valeur optimale primale: ', primal);
   % fprintf('La valeur optimale duale: %.10f\n', dual);
   % fprintf('Le nombre d'itération externe est', kext);

 
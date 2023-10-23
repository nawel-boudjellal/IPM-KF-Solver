function [s]=pointinit(A,b1)    
%    Le point initiale primale vérifie Ax=b1,x>=0
% A=[1 1 1 1;1 1 0 -3];
% b1=[1;0.5];
% c=[1;2;3;4];
n=size(A,2);
% donner la point initial x0>0
x=ones(n,1);
x(n+1)=1;
%%%%
An=[A,b1-A*x(1:n)];
c1=[zeros(n,1);1];
n=n+1;
g=0;
% la procedure de Karmarkar généralisé
 epsilon=10^(-4);
 %alpha=1/4;
 alpha=n/(2*n-1);
while(c1'*x-g>epsilon)
    Dk=diag(x);
    Ak=An*Dk;
     Bk=[Ak,-b1;ones(1,n+1)];
     Dkc=[Dk*c1;-g];
      z=Bk*Dkc;
      ZZ=Bk*(Bk)';
     R=chol(ZZ);
     %CC=det(ZZ)
      u=gauss1(R',z);
     T=gauss2(R,u);
    pk=Dkc-Bk'*T;
    dk=pk/norm(pk,2);
     r=1/(sqrt(n*(n+1)));
    y=ones(n+1,1)-alpha*r*dk;
    x=Dk*y(1:n)/y(n+1);
 end

s=x(1:(n-1));
% c1'*x
% An*x-b1
% A*x(1:(n-1))
end
% Central Trajectory Method via a new kernel function for convex (Linear) optimization
%clc
%clear all
tic
theta=0.99;
beta=0.99;
 epsilon=10^(-4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Netlib example
load("lp_scagr7.mat")
A = Problem.A; 
a=full(A);
b = Problem.b;
b=full(b);
c = Problem.aux.c;
c=full(c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(A);
%The initial point is found by executing PrimalDualPL.m file
% x=zeros(n,1);
% z=zeros(n,1);
 % y=zeros(m,1);
 p=log(n)/2;
tau=sqrt(n);
 mu=x'*z/n;
kext=0;
kint=0;
 psi=inline('(t^2-1)/2+log(2/(1+t))-(exp(1/(1+t)-0.5)/(2*p))-(1/(2*p))');
psiv=zeros(n,1); 
psi1=inline('t-(1/(1+t))-(2*(exp(1/(1+t)-0.5)/(2*p))/((1+t)^2))');
psi1v=zeros(n,1);
 while x'*z>epsilon
     %external iteration
     mu=(1-theta)*mu;
     v=sqrt(x.*z/mu);
     for i=1:n
         psiv(i)=psi(p,v(i));
     end
     phi=sum(psiv);
   
     while phi>tau
       %Interne iteration
       for i=1:n
           psi1v(i)=psi1(p,v(i));
       end
       r=-psi1v;
       s=zeros(m,1);
       t=zeros(n,1);
       a1=1/mu*a*(diag(v)^(-1))*diag(x);
       M=a1*a1';
       R=chol(M);
       s1=s-a1*(r+t);
       u=R'\s1;
       deltay=R\u;
       %
       dx=a1'*deltay+(r+t);
       dz=-a1'*deltay-t;
       deltax=x.*dx./v;
       deltaz=z.*dz./v;
       delta=1/2*norm(psi1v,2);
      alphax=depl(x,deltax);
      alphaz=depl(z,deltaz);
     alpha=beta*min(alphax,alphaz);
       x=x+alpha*deltax;
       y=y+alpha*deltay;
       z=z+alpha*deltaz;
       kint=kint+1;
       v=sqrt(x.*z/mu);
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
% fprintf('Gap value is', saut);
% fprintf('optimal solution is', x,y,z)
  %fprintf('The primal optimal value: ', primal);
   % fprintf('The dual optimal value: %.10f\n', dual);
   % fprintf('The number of external iterations is', kext);
 % fprintf('The number of interne iterations is', kint);
 

% Initial solution strictly feasible primal dual of linear problem
% (LP)
           %Ax=b,x>=0
           %A'y+z=c,z>=0
clc
clear all
%tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Example Netlib
%Name of problem lp_scagr7.mat
load("lp_scagr7.mat")
A = Problem.A; 
A=full(A);
b = Problem.b;
b=full(b);
c = Problem.aux.c;
c=full(c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(A);
B=[A,zeros(m,m),zeros(m,m),zeros(m,n);zeros(n,n),A',-A',eye(n)];
S=[b;c];
w=pointinit(B,S);
%%%%%%
x=w(1:n);
y=w(n+1:m+n)-w(m+n+1:2*m+n);
z=w(2*m+n+1:2*(m+n));
%input('le point initial')
disp(x)
disp(y)
disp(z)
%toc

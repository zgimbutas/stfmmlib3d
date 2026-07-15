function [x,ier]=gmres_simple(A,b,tol,maxit,x0)
%GMRES algorithm.
%
%  This subroutine solves a complex linear system Ax=b by means
%  of GMRES algorithm.
%
%  X = GMRES_SIMPLE(A,b);
%  X = GMRES_SIMPLE(A,b,tol);
%  X = GMRES_SIMPLE(A,b,tol,maxit);
%  X = GMRES_SIMPLE(A,b,tol,maxit,x0);
%
%  Input parameters:
%
%  A - the matrix of the system (or the function handle to evaluate A(x) )
%  b - the right hand side 
%  tol - the required accuracy
%  maxit - the maximum number of iteration permitted
%  x0 - the initial guess
%
%  Output parameters:
%
%  x - the solution of the system
%  ier - error return code
%     ier=0 normal execution of the subroutine
%     ier=4 means that the maximum number iterations maxit
%           has been reached without achieving the required accuracy tol
%     ier=8 means that the errors failed to decrease before the maximum
%           number of iterations has been reached or the required accuracy
%           eps has been reached
%

[n,m]=size(b);

if( nargin < 3 ), tol = 1e-6; end
if( nargin < 4 ), maxit = min(n,20); end

Ae=zeros(n,maxit);
e=zeros(n,maxit);

ier=0;

% A is a matrix, construct function handle 
if( isnumeric(A) ), A = (@(x) A*x); end

if( nargin == 5),
x=x0;
r=b-A(x);
else
x=zeros(n,1);
r=b;
end
er0=norm(r,2);

% Randomize the initial solution vector
%x=rand(n,1);
%r=b-A(x);
%er0=norm(r,2);

for i=1:maxit

% recalculate residual every 10 steps
if( mod(i,10) == 0), r=b-A(x); end;

% initial new direction guess
w=A(r);

% orthogonalize the new direction to the preceding ones
wp=w;
rp=r;
for l=1:1
for j=1:i-1
  d=Ae(:,j)'*wp;
  wp=wp-d*Ae(:,j);
  rp=rp-d*e(:,j);
end
end

% normalize and store the new direction
d=1/sqrt(wp'*wp);
Ae(:,i)=d*wp;
e(:,i)=d*rp;

% reorthogonalize the new direction one more time
for l=1:1
for j=1:i-1
  d=Ae(:,j)'*Ae(:,i);
  Ae(:,i)=Ae(:,i)-d*Ae(:,j);
  e(:,i)=e(:,i)-d*e(:,j);
end
end

% normalize and store the new direction
d=1/sqrt(Ae(:,i)'*Ae(:,i));
Ae(:,i)=Ae(:,i)*d;
e(:,i)=e(:,i)*d;

% update the residual and solution
d=Ae(:,i)'*r;
r=r-d*Ae(:,i);

% the residual stopped decreasing, abort
er1=norm(r,2);
if( er1 >= er0 ),
  ier=8; return; 
end;
er0=er1;

x=x+d*e(:,i);

%%%norm(r,2)
%%%fprintf('iter: %d, norm(r,2): %f\n',i,norm(r,2));
%%%fprintf('iter: %d, rms=norm(r,2)/sqrt(n): %17.12f\n',i,norm(r,2)/sqrt(n));
fprintf('iter: %d, rel=norm(r,2)/norm(b,2): %17.12e\n',i,norm(r,2)/norm(b,2));

%%%if( i == 1 ), norm1 = norm(r,2); end;
if( i == 1 ), norm1 = norm(b,2); end;
if( norm(r,2) < tol*norm1 ), return; end;

end

% no convergence, abort
ier=4;

function [x,ier]=bicgstab_simple(A,b,tol,maxit,x0)
% BiCG(stab) algorithm
%
%  This subroutine solves a complex linear system Ax=b by means
%  of BiCG(stab) algorithm. It is equivalent to bicgstab_l with ell=1.
%
%  X = BICGSTAB_SIMPLE(A,b);
%  X = BICGSTAB_SIMPLE(A,b,tol);
%  X = BICGSTAB_SIMPLE(A,b,tol,maxit);
%  X = BICGSTAB_SIMPLE(A,b,tol,maxit,x0);
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
%     ier=8 breakdown (a denominator vanished; no look-ahead recovery)
%
%  Cost: 2 matrix-vector products per iteration.
%
%  Reference: H.A. Van der Vorst, "Bi-CGSTAB: A fast and smoothly
%  converging variant of Bi-CG for the solution of nonsymmetric linear
%  systems", SIAM J. Sci. Stat. Comput., 13(2):631-644, 1992.

[n,m]=size(b);

if( nargin < 3 ) tol = 1e-6; end
if( nargin < 4 ) maxit = min(n,20); end

ier=0;

% A is a matrix, construct function handle 
if( isnumeric(A) ), A = (@(x) A*x); end

% initialize 
if( nargin == 5 ),
xk=x0;
else
xk=zeros(n,1);
end
axk=A(xk);

if( n == 1 )
  x = b / A(1);
  if( ~isfinite(x) ), x=xk; ier=8; end   % singular scalar operator
  return;
end


% find the first direction
er = 0;
er1 = 0;

ek=b-axk;
ekm1 = ek;

bek = ek;
bekm1 = bek;

norm1=norm(b,2);

% zero RHS: the solution is zero
if( norm1 == 0 ), x=zeros(n,1); return; end

% initial guess already converged: skip recurrence (avoids 0/0 NaN)
if( norm(ek,2) <= tol*norm1 ), x=xk; return; end

er1=(bek'*ek);



% bicgstab algorithm
for i=1:maxit

aekm1=A(ekm1);

cd1=bek'*ek;
cd2=bek'*aekm1;

alpha=cd1/cd2;
if( ~isfinite(alpha) ), x=xk; ier=8; return; end

esk=ek-alpha*aekm1;

% intermediate residual check: if s=r-alpha*v is already small the
% alpha half-step alone converged; stop here (also avoids 0/0 in omega
% when esk=0, e.g. the exact solution is reached in one step)
if( norm(esk,2) < tol*norm1 ), x=xk+alpha*ekm1; return; end

etk=A(esk);

cd1=etk'*esk;
cd2=etk'*etk;

omega=cd1/cd2;
if( ~isfinite(omega) ), x=xk; ier=8; return; end

% update the solution
xk=xk+alpha*ekm1+omega*esk;

% update the residual
ek=esk-omega*etk;

er=(ek'*ek);
er2=(bek'*ek);

%%%norm(r,2)
%%%fprintf('iter: %d, norm(r,2): %f\n',i,norm(r,2));
%%%fprintf('iter: %d, rms=norm(r,2)/sqrt(n): %17.12f\n',i,er/sqrt(n));
fprintf('iter: %d, rel=norm(r,2)/norm(b,2): %17.12e\n',i,sqrt(er)/norm(b,2));

% check convergence before forming the next direction, so a converged
% step is reported even if the next beta would break down
if( sqrt(er) < tol*norm1 ), x=xk; return; end;

% prepare the next search direction
beta=(er2/er1)*(alpha/omega);
er1=er2;
if( ~isfinite(beta) ), x=xk; ier=8; return; end

ekm1=ek+beta*(ekm1-omega*aekm1);

end

% % no convergence, abort, return the solution
x=xk;
ier=4;

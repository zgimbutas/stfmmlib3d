function x=bicgstab_simple(A,b,tol,maxit)
% BiCG(stab) algorithm


[n,m]=size(b);

if( nargin < 3 ) maxit = min(n,20); end
if( nargin < 4 ) tol = 1e-6; end

% initialize 
xk=zeros(n,1);
axk=A(xk);

if( n == 1 )
  xk = b/axk;
  axk = b;
end


% find the first direction
er = 0;
er1 = 0;

ek=b-axk;
ekm1 = ek;

bek = ek;
bekm1 = bek;

er0=sqrt(ek'*ek);
er1=sqrt(bek'*ek);



% bicgstab algorithm
for i=2:maxit

aekm1=A(ekm1);

cd1=bek'*ek;
cd2=bek'*aekm1;

alpha=cd1/cd2;

esk=ek-alpha*aekm1;
etk=A(esk);

cd1=etk'*esk;
cd2=etk'*etk;

omega=cd1/cd2;

% update the solution
xk=xk+alpha*ekm1+omega*esk;

% update the residual
ek=esk-omega*etk;

er=sqrt(ek'*ek);
er2=sqrt(bek'*ek);

beta=(er2/er1)*(alpha/omega);
er1=er2;

ekm1=ek+beta*(ekm1-omega*aekm1);

fprintf('iter: %d, norm(r,2)/sqrt(n): %17.12f\n',i,er/sqrt(n));
fprintf('iter: %d, norm(r,2)/norm(b,2): %17.12f\n',i,er/norm(b,2));

%%%if( i == 1 ) er0=er; end;
if( i == 1 ) er0=norm(b,2); end;
if( er < tol*er0 ), break; end;

end

% return the solution
x=xk;

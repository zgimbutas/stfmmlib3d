%
%  Test Stokes particle FMMs in R^3
%

addpath('/home/zng2/todo/fmmlib3d-1.2/matlab');

nsources = [5000 5000 50000 500000]
iprecs = [-1 0 1 2]

for i=1:size(nsources,2)
    
nsource = nsources(i);

source = zeros(3,nsource);

idist=3;

if( idist == 1 ),
theta=rand(1,nsource)*pi;
phi=rand(1,nsource)*2*pi;
source(1,:)=.5*cos(phi).*sin(theta);
source(2,:)=.5*sin(phi).*sin(theta);
source(3,:)=.5*cos(theta) + 2;
end

if( idist == 2 ),
source(1,:)=rand(1,nsource);
source(2,:)=rand(1,nsource);
source(3,:)=rand(1,nsource) + 2;
end

if( idist == 3 ),
phi=rand(1,nsource)*2*pi;
source(1,:)=.5*cos(phi);
source(2,:)=.5*sin(phi);
source(3,:)=rand(1,nsource) + 2;
end

%scatter3(source(1,:),source(2,:),source(3,:))
%plot3(source(1,:),source(2,:),source(3,:))

for j=1:size(iprecs,2)

iprec = iprecs(j);

%
%  Timings
%


ifsingle=0;
sigma_sl = rand(3,nsource);
ifdouble=1;
sigma_dl = rand(3,nsource);
sigma_dv = rand(3,nsource);

ifpot = 1;
ifgrad = 0;

'Stokes particle FMM in R^3';

tic
[U]=stfmm3dpart(iprec,nsource,source,ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,ifpot,ifgrad);
total_time=toc;
total_time_fmm=total_time;
speed=nsource/total_time;

'Stokes particle direct evaluation in R^3';

k = 100;
if( nsource < 10000), k = 1; end;
if( nsource > 100000), k = 1000; end;

ms = nsource / k;
tic
[F]=st3dpartdirecttime(ms,nsource,source,ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,ifpot,ifgrad);
total_time=toc;
total_time_direct=total_time*k;
speed=nsource/total_time;


if( ifpot ), U.pot=U.pot/(4*pi); end
if( ifgrad ), U.grad=U.grad/(4*pi); end

if( ifpot ), F.pot=F.pot/(4*pi); end
if( ifgrad ), F.grad=F.grad/(4*pi); end



if( ifpot ), U.pot=U.pot(:,1:ms); end
if( ifgrad ), U.grad=U.grad(:,1:ms); end

if( ifpot ), F.pot=F.pot(:,1:ms); end
if( ifgrad ), F.grad=F.grad(:,1:ms); end


if( ifpot ),
rel_error_pot = norm((U.pot - F.pot),2)/norm((F.pot),2);
end

if( ifgrad ),
rel_error_grad = norm(U.grad - F.grad,2)/norm(F.grad,2);
end

if( ifpot ),
  fprintf('%5d & %2d & %5.3f & %5.3f & %6.3e \\\\ \n', nsource, iprec, total_time_fmm, total_time_direct, rel_error_pot );
end

if( ifgrad ),
  fprintf('%5d & %2d & %5.3f & %5.3f & %6.3e \\\\ \n', nsource, iprec, total_time_fmm, total_time_direct, rel_error_grad );
end

end

end

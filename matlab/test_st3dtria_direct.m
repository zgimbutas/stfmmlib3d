%
%  Test Stokes triangle FMMs in R^3
%

%
%  Retrieve flat triangulation
%

geom_type = 2;
filename_geo = 'sphere180.a.tri';
filename_geo = 'sphere720.a.tri';
%filename_geo = 'sphere2880.a.tri';
%filename_geo = 'sphere11520.a.tri';
%filename_geo = 'sphere20480.a.tri';

fid = fopen(filename_geo,'r');

nverts=0;
nfaces=0;
[nverts] = fscanf(fid,'%d',1);
[nfaces] = fscanf(fid,'%d',1);

verts=zeros(3,nverts);
ifaces=zeros(3,nfaces);

[verts] = fscanf(fid,'%f',[3,nverts]);
[ifaces] = fscanf(fid,'%d',[3,nfaces]);

fclose(fid);

%
%  refined rectangle
%
%test2a
%filename_geo='rectangle'

nverts,nfaces

%
%  create triangle vertex and normal arrays
%

ntri = nfaces;
triangles = zeros(3,3,ntri);

for i=1:ntri
	
%triangles(1:3,1,i) = verts(1:3,ifaces(1,i));
%triangles(1:3,2,i) = verts(1:3,ifaces(2,i));
%triangles(1:3,3,i) = verts(1:3,ifaces(3,i));

triangles(1:3,1:3,i) = verts(1:3,ifaces(1:3,i));

end

%
%  build triangle normals and area vector
%
trianorm = zeros(3,ntri);
triaarea = zeros(1,ntri);

for i=1:ntri

vec1 = triangles(1:3,2,i) - triangles(1:3,1,i);
vec2 = triangles(1:3,3,i) - triangles(1:3,1,i);

trianorm(1:3,i) = cross(vec1,vec2);
triaarea(i) = norm(trianorm(1:3,i))/2;

trianorm(1:3,i) = trianorm(1:3,i)/norm(trianorm(1:3,i));

end

%%% sum(triaarea), 4*pi

%
%  centroids
%

source = sum(triangles,2)/3;
source = reshape(source,3,ntri);

%plot3(source(1,:),source(2,:),source(3,:))

nsource = ntri

%
%  timings
%

ifsingle=1;
sigma_sl = rand(3,ntri);
ifdouble=0;
sigma_dl = rand(3,ntri);
dipvec = trianorm;

%%sigma_dl = zeros(3,ntri);
%%sigma_dl(2,:) = 1;
%%sigma_dl(2,:) = 1*cos(.1*source(3,:));
%%sigma_dl = cross(sigma_dl,cross(sigma_dl,trianorm));

ifsingle
ifdouble
ifpot = 1
ifgrad = 0

ntarget = ntri
target = source(:,1:ntri);
target(1,:) = target(1,:) + 10;

%%h=1e-4
%%target = [source + h*trianorm source - h*trianorm];
%%target = [source + h*trianorm];
%%target = [source - h*trianorm];

[ndim,ntarget] = size(target)
ifpottarg = 1
ifgradtarg = 0
ntarget
%ntarget = 0


disp('')
'Stokes triangle target FMM in R^3'

tic
iprec=0
[U]=stfmm3dtria(iprec,ntri,triangles,trianorm,source,ifsingle,sigma_sl,ifdouble,sigma_dl,ifpot,ifgrad,ntarget,target,ifpottarg,ifgradtarg);
total_time=toc


'Stokes triangle target direct evaluation in R^3'
tic
[F]=st3dtriadirect(ntri,triangles,trianorm,source,ifsingle,sigma_sl,ifdouble,sigma_dl,ifpot,ifgrad,ntarget,target,ifpottarg,ifgradtarg);
total_time=toc


if( ifpot ), U.pot=U.pot/(4*pi); end
if( ifpot ), U.pre=U.pre/(4*pi); end
if( ifgrad ), U.grad=U.grad/(4*pi); end

if( ifpot ), F.pot=F.pot/(4*pi); end
if( ifpot ), F.pre=F.pre/(4*pi); end
if( ifgrad ), F.grad=F.grad/(4*pi); end

if( ifpot ),
%rms_pot = norm((F.pot),2)/sqrt(nsource)
rms_error_pot = norm((U.pot - F.pot),2)/sqrt(nsource)
rms_error_pre = norm((U.pre - F.pre),2)/sqrt(nsource)
rel_error_pot = norm((U.pot - F.pot),2)/norm((F.pot),2)
rel_error_pre = norm((U.pre - F.pre),2)/norm((F.pre),2)
end

if( ifgrad ),
%rms_grad = norm(reshape(F.grad,9,nsource),2)/sqrt(nsource)
rms_error_grad = norm(reshape(U.grad - F.grad,9,nsource),2)/sqrt(nsource)
rel_error_grad = norm(reshape(U.grad - F.grad,9,nsource),2)/norm(reshape(F.grad,9,nsource),2)
end
%%%break;

if( ifpottarg ), U.pottarg=U.pottarg/(4*pi); end
if( ifpottarg ), U.pretarg=U.pretarg/(4*pi); end
if( ifgradtarg ), U.gradtarg=U.gradtarg/(4*pi); end

if( ifpottarg ), F.pottarg=F.pottarg/(4*pi); end
if( ifpottarg ), F.pretarg=F.pretarg/(4*pi); end
if( ifgradtarg ), F.gradtarg=F.gradtarg/(4*pi); end

if( ifpottarg ),
%rms_pottarg = norm((F.pottarg),2)/sqrt(nsource)
rms_error_pottarg = norm((U.pottarg - F.pottarg),2)/sqrt(ntarget)
rms_error_pretarg = norm((U.pretarg - F.pretarg),2)/sqrt(ntarget)
norm_pottarg = norm((F.pottarg),2)
rel_error_pottarg = norm((U.pottarg - F.pottarg),2)/norm((F.pottarg),2)
rel_error_pretarg = norm((U.pretarg - F.pretarg),2)/norm((F.pretarg),2)
end

if( ifgradtarg ),
rms_gradtarg = norm(reshape(F.gradtarg,9,ntarget),2)/sqrt(ntarget)
rms_error_gradtarg = ...
    norm(reshape(U.gradtarg - F.gradtarg,9,ntarget),2)/sqrt(ntarget)
rel_error_gradtarg = ...
    norm(reshape(U.gradtarg - F.gradtarg,9,ntarget),2)/ ...
    norm(reshape(F.gradtarg,9,ntarget),2)
end
%%%break;


disp('')
'Stokes triangle FMM in R^3'

ifpottarg = 0
ifgradtarg = 0


tic
iprec=0
[U]=stfmm3dtria(iprec,ntri,triangles,trianorm,source,ifsingle,sigma_sl,ifdouble,sigma_dl,ifpot,ifgrad);
total_time=toc

'Stokes triangle direct evaluation in R^3'

tic
[F]=st3dtriadirect(ntri,triangles,trianorm,source,ifsingle,sigma_sl,ifdouble,sigma_dl,ifpot,ifgrad);
total_time=toc



if( ifpot ), U.pot=U.pot/(4*pi); end
if( ifgrad ), U.grad=U.grad/(4*pi); end

if( ifpot ), F.pot=F.pot/(4*pi); end
if( ifgrad ), F.grad=F.grad/(4*pi); end

if( ifpot ),
%rms_pot = norm((F.pot),2)/sqrt(nsource)
rms_error_pot = norm((U.pot - F.pot),2)/sqrt(nsource)
end

if( ifgrad ),
%rms_grad = norm(reshape(F.grad,9,nsource),2)/sqrt(nsource)
rms_error_grad = norm(reshape(U.grad - F.grad,9,nsource),2)/sqrt(nsource)
end
%%%break;


%filename_out=['output_' filename_geo '.mat']
%save(filename_out,'-v6')

%plot_solution
%plot_displacement
%plot_slip



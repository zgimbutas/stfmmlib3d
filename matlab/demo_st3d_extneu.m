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
%filename_geo = 'cube_r1.a.tri';
%filename_geo = 'cube_r2.a.tri';
%filename_geo = 'cube_r3.a.tri';
%filename_geo = 'cube_r4.a.tri';
%filename_geo = '../../Data/cube768.a.tri';
%filename_geo = '../../Data/cube3072.a.tri';

[verts,ifaces,nverts,nfaces] = atriread(filename_geo);
nverts,nfaces


%
%  Construct triangle vertex, normal, area, and centroid arrays
%

ntri = nfaces;
[triangles,trianorm,triaarea,source]=atriproc(verts,ifaces);

%plot3(source(1,:),source(2,:),source(3,:),'*')

nsource = ntri;

%
%  Initialize matrix multiplication routine
%


nsource = ntri;

A.alpha = 1/(2*pi);
A.iprec = 0;

A.ntri = ntri;
A.triangles = triangles;
A.trianorm = trianorm;
A.source = source;

smu = 1
A.smu = smu;


%
%  Construct the test right hand side
%
%%%rhs = ones(3,nsource);

ntest_source = 1;
iftest_charge = 1;
test_source = [.1;.2;.3];
test_charge = [0;0;1];
iftest_dipole = 0;
test_dipstr = [0;1;0];
test_dipvec = [1;0;0];

ifpot = 0;
ifgrad = 0;
ifpottarg = 1;
ifgradtarg = 1;

[U]=st3dpartdirect(ntest_source,test_source,iftest_charge,test_charge,iftest_dipole,test_dipstr,test_dipvec,ifpot,ifgrad,nsource,source,ifpottarg,ifgradtarg);

%%%rhs = U.pottarg;
rhs = st3dtraction(U.pretarg,U.gradtarg,A.trianorm);


M = sqrt(triaarea)';
%
%  Call the solver
%
'Exterior Neumann solver for the Stokes equation in R^3'
tic
rhs0 = reshape(rhs,1,3*nsource);
sol0 = gmres_simple(@(x) st3dmultfmmflat_neu(A,x), rhs0.', 1e-3, 20);
sol0 = sol0.';
sol = reshape(sol0,3,nsource);
time_gmres=toc


%
%  Evaluate potential and field directly
%
ntarget = 1;
target = [10;-20;30];


ifpottarg = 1;
ifgradtarg = 0;
[U]=st3dpartdirect(ntest_source,test_source,iftest_charge,test_charge,iftest_dipole,test_dipstr,test_dipvec,ifpot,ifgrad,ntarget,target,ifpottarg,ifgradtarg);

pot = U.pottarg;


%
%  Evaluate potential and field via the solution of integral equation
%
tic
iprec=0;

charge = sol;
dipstr = zeros(3,nsource);
dipvec = trianorm;

ifcharge=1;
ifdipole=0;
ifpot=1;
ifgrad=0;
ifpottarg=0;
ifgradtarg=0;
[P]=stfmm3dtria(iprec,nsource,triangles,trianorm,source,...
    ifcharge,charge,ifdipole,dipstr,ifpot,ifgrad,...
    ntarget,target,ifpottarg,ifgradtarg);

charge = sol;
dipstr = A.alpha * P.pot;
dipvec = trianorm;

ifcharge=1;
ifdipole=1;
ifpot=0;
ifgrad=0;
ifpottarg=1;
ifgradtarg=0;
[U]=stfmm3dtria(iprec,nsource,triangles,trianorm,source,...
    ifcharge,charge,ifdipole,dipstr,ifpot,ifgrad,...
    ntarget,target,ifpottarg,ifgradtarg);

time_postproc=toc

%
%  Finally, print the results
%
'Potential and field at the targets'
pot,U.pottarg 

rel_error_pot=norm((pot-U.pottarg)./pot)


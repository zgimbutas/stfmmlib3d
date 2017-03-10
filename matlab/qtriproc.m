function [triaquad,trianorm,triaarea,centroid,triatang1,triatang2]=qtriproc(verts,ifaces)
%QTRIPROC Process triangulations in Cart3d format (quadratic).
%
%  [triaquad]=qtriproc(verts,ifaces);
%
%  [triaquad,trianorm,triaarea,centroid]=qtriproc(verts,ifaces);
%
%  [triaquad,trianorm,triaarea,centroid,triatang1,triatang2]=
%           qtriproc(verts,ifaces);
%
%  Input parameters:
%
%  verts - real(3,nverts): array of triangulation vertices
%  ifaces - integer(6,nfaces): indices of triangle vertices
%
%  Output parameters:
%
%  triaquad - real(3,6,ntri): array of triangle vertex and midpoint coordinates 
%  trianorm - real(3,nsource): triangle normals at centroids
%  triaarea - real(nsource): triangle area elements at centroids
%  centroid - real(3,nsource): triangle centroids
%  triatang1 - real(3,nsource): triangle tangents at centroids (first set)
%  triatang2 - real(3,nsource): triangle tangents at centroids (second set)
%
%  Note: the first set of tangent vectors is (\partial xyz/\partial u).
%

%
%  Construct triangle vertex array
%
nverts=size(verts,2);
nfaces=size(ifaces,2);

ntri = nfaces;
triaquad = zeros(3,6,ntri);

%for i=1:ntri	
%  triaquad(1:3,1:6,i) = verts(1:3,ifaces(1:6,i));
%end

triaquad(1:3,1:6,1:ntri) = reshape(verts(1:3,ifaces(1:6,1:ntri)),3,6,ntri);


if( nargout > 1 ),
%
%  Parametrization constants
%
%       ... setup a quadratic triangle in R^3
%
%
%            2     
%          .   . 
%         C     B
%        .       .
%       0 .. A .. 1
%
%
x0=squeeze(triaquad(1:3,1,:));
x1=squeeze(triaquad(1:3,2,:));
x2=squeeze(triaquad(1:3,3,:));
xa=squeeze(triaquad(1:3,4,:));
xb=squeeze(triaquad(1:3,5,:));
xc=squeeze(triaquad(1:3,6,:));

xu=-3*x0+4*xa-x1;
xv=-3*x0+4*xc-x2;
xuu=x0-2*xa+x1;
xuv=x0-xa+xb-xc;
xvv=x0-2*xc+x2;

%
%  Triangle centroids
%
u=1/3;
v=1/3;

centroid=x0+u*xu+v*xv+2*(u*(u*xuu+v*xuv)+v*(u*xuv+v*xvv));

%
%  Construct triangle normals
%
trianorm = zeros(3,ntri);
triaarea = zeros(1,ntri);

vec1 = xu+4*(u*xuu+v*xuv);
vec2 = xv+4*(u*xuv+v*xvv);

trianorm = cross(vec1,vec2);
ds = sqrt(sum(trianorm.^2,1));

trianorm = trianorm ./ repmat(ds,3,1);

%  Triangle area element at the centroid
triaarea = ds/2;
end


if( nargout > 4 ),
%
%  Construct tangent vectors
%
dt = sqrt(sum(vec1.^2,1));
triatang1 = vec1 ./ repmat(dt,3,1);
triatang2 = cross(trianorm,triatang1);
end

function [verts,ifaces,nverts,nfaces]=qtriread(filename)
%QTRIREAD Retrieve triangulations in Cart3d from a file.  (quadratic)
%
%  [verts,ifaces,nverts,nfaces]=qtriread(filename);
%
%  Input parameters:
%
%  filename - input file name.
%
%  Output parameters:
%
%  verts - real(3,nverts): array of triangulation vertices
%  ifaces - integer(6,nfaces): indices of triangle vertices
%

%
%  Retrieve quadratic triangulation
%

geom_type = 3;

fid = fopen(filename,'r');

nverts=0;
nfaces=0;
[nverts] = fscanf(fid,'%d',1);
[nfaces] = fscanf(fid,'%d',1);

verts=zeros(3,nverts);
ifaces=zeros(3,nfaces);

[verts] = fscanf(fid,'%f',[3,nverts]);
[ifaces] = fscanf(fid,'%d',[6,nfaces]);

fclose(fid);

%%%nverts,nfaces


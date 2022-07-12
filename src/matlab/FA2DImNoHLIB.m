% input:
%---X,Y: input array
%---[rsc rsc]: array with two elements
function [qsavex,qsavey] = FA2DImNoHLIB(X,Y,wt,niter,mode,lscale)

% Get the dimension sizes
sx = size(X,1);
sy = size(X,2);

% Get the nearest lesser power of 2
sx2 = 2*pow2(floor(log2(sx)));
sy2 = 2*pow2(floor(log2(sy)));

% resize input and target images to the input size
xx = myimresize(X,[sx2 sy2]); 
yy = myimresize(Y,[sx2 sy2]);

[qtx,qty] = FA2DC4V2_simple(xx,yy,niter,wt,lscale,mode);
%H=ones(size(xx));
%[qtx,qty,qe1,qe2] = FA2DC4V2LIB(xx,yy,H,[],[],[],1.0,niter,wt,lscale,mode);

qsavex = myimresize(qtx,[size(X,1) size(X,2)])*size(X,2)/sx2;
qsavey = myimresize(qty,[size(X,1) size(X,2)])*size(X,1)/sy2;

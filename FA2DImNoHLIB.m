% input:
%---X,Y: input array
%---[rsc rsc]: array with two elements
function [qsavex,qsavey] = FA2DImNoHLIB(X, Y,rsc,wt,niter,mode,lscale)
xx = myimresize(X,[rsc rsc]); 
yy = myimresize(Y,[rsc rsc]);
H=ones(rsc);
[qtx,qty,qq,err]=FA2DC4V2LIB(xx,yy,H,[],[],[],1,1,...
                           niter,0,wt,lscale,mode,1);


qsavex = myimresize(qtx,[size(X,1) size(X,2)])*size(X,2)/rsc;
qsavey = myimresize(qty,[size(X,1) size(X,2)])*size(X,1)/rsc;

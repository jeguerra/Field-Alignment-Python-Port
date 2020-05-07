function [qsavex,qsavey] = FA2DImNoH(X, Y,rsc,wt,niter,mode,lscale)
xx = imresize(X,[rsc rsc]); yy = imresize(Y,[rsc rsc]);
[qtx,qty,qq,err]=FA2DC4V2(xx,yy,ones(rsc),[],[],[],1,1,...
                           niter,0,wt,lscale,mode,1);

qsavex = imresize(qtx,[size(X,1) size(X,2)])*size(X,2)/rsc;
qsavey = imresize(qty,[size(X,1) size(X,2)])*size(X,1)/rsc;


  


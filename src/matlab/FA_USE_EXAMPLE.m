%close all;
clear all;
addpath(genpath('/scratch1/portfolios/BMC/fast/IJankov/FASTv4.2a/FA2D'));

for count=1:120;
  X=eval(['readgfsbin(''/scratch1/portfolios/BMC/fast/IJankov/FASTv4.2a/experiment/test' num2str(count) '/fldA.bin'')']); %forecast
  Y=eval(['readgfsbin(''/scratch1/portfolios/BMC/fast/IJankov/FASTv4.2a/experiment/test' num2str(count) '/fldB.bin'')']); %analysis
  mX = max([X(:); Y(:)]); 
  iter = 2^10;
  lscale = 0;
  qtx = zeros(size(X)); qty = qtx;
  [qx,qy] = FA2DImNoHLIB(X,Y,128,.1,iter,2,lscale);
  X = X.*mX;Y = Y.*mX;
  FX = advect(X,qx,qy);

  Xn(:,:,count)=X; clear X
  Yn(:,:,count)=Y; clear Y
  FXn(:,:,count)=FX; clear FX
  
  qtxn(:,:,count)=qtx; clear qtx
  qtyn(:,:,count)=qty; clear qty
  qxn(:,:,count)=qx; clear qx  
  qyn(:,:,count)=qy; clear qy
  
  save counter.txt count -ascii
end

Xmean=mean(Xn,3);
Ymean=mean(Yn,3);
FXmean=mean(FXn,3);

qtxmean=mean(qtxn,3);
qtymean=mean(qtyn,3);
qxmean=mean(qxn,3);
qymean=mean(qyn,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=Xmean;
Y=Ymean;
FX=FXmean;
qtx=qtxmean;
qty=qtymean;
qx=qxmean;
qy=qymean;

save DATA_4plots X Y FX qtx qty qx qy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=figure;
subplot(221);
imagesc(Y-X,[-500 500]); colorbar;
axis('tight'); axis('xy');
title('Original Amplitude Error','FontSize',12);

subplot(222);sy = size(qx,1); sx = size(qx,2);
%quiver(qx(1:5:sy,1:5:sx),...
%    qy(1:5:sy,1:5:sx),1.2,'k');
quiver(qx(1:5:sy,1:5:sx),...
    qy(1:5:sy,1:5:sx),2.5,'k');
axis('tight');axis('xy');axis('off');
title('Estimated Displacement Component','FontSize',12);


subplot(223);
imagesc(Y-FX,[-500 500]); colorbar;
axis('tight');
axis('xy');axis('off');
title('Residual Amplitude Error','FontSize',12);

subplot(224);imagesc(FX-X,[-500 500]); colorbar;
axis('tight');axis('xy');axis('off');
title('Original Amplitude Error Explained by FA','FontSize',12);





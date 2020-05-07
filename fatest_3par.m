%close all;
clear all;
t = cputime;
%Field alignment
%input:
%---two fields---
%A and B
%output:
%--C--> qx, qy : vector field
%--D align B to A
%---
%addpath(genpath('/scratch1/portfolios/BMC/fast/IJankov/FASTv4.2a/FA2D'));

%for count=1:120;
%input files------------------
%dir = ...    filename = ...
    fileID1 = fopen('/scratch3/BMC/det/Jie.Feng/Jing/readoutmem/GH_20mem_500hPa_0to360andS90toN90_120h.dat')
%output files--------------------
%-----------------------------------
    precision = 'float';
    nlon = 360;
    nlat = 181;
    hmem0 = fread(fileID1,nlon*nlat*2,precision);   % read the file Yall(64*32*2)
    hmem = reshape(hmem0,[nlon,nlat,2]);      
%    sz = size(Yall1);
%--------------------------
    hmemT = zeros(nlat,nlon,2,'single');
    hmemT(:,:,1) = transpose(hmem(:,:,1));  %matrix transposed before input to subroutine
    hmemT(:,:,2) = transpose(hmem(:,:,2));  %matrix transposed before input to subroutine
   %============================== 
   fprintf('Start Calculation');
    %mX = max([X(:); Y(:)]); 
    iter = 2^10;
    lscale = 1;   %change from 1~10,   1 is the best 
%=============================
    qxyT = zeros(nlat,nlon,2,'single');
    qxy = zeros(nlon,nlat,2,'single');
%------use the function-------
%------spd1T: B------
%------spd0T: A------
%2->1
    %[qxyT(:,:,1),qxyT(:,:,2)] = FA2DImNoHLIB(hmemT(:,:,2),hmemT(:,:,1),128,.1,iter,2,lscale);
     [qxyT(:,:,1),qxyT(:,:,2)] = FA2DImNoHLIB(hmemT(:,:,1),hmemT(:,:,2),128,.1,iter,2,lscale);
	fclose(fileID1);
    %fileID3 = fopen('/scratch3/BMC/det/Jie.Feng/QGmodel/data/Errest/gridmean/est/par3/afterFA/displacement_qx.dat','w')
    %fileID4 = fopen('/scratch3/BMC/det/Jie.Feng/QGmodel/data/Errest/gridmean/est/par3/afterFA/displacement_qy.dat','w')
	qxy(:,:,1) = transpose(qxyT(:,:,1));    %output after transopose
	qxy(:,:,2) = transpose(qxyT(:,:,2));
%	fwrite(fileID3,qx,precision)  %u wind
%	fwrite(fileID4,qy,precision)  %v wind
    %X = X.*mX;Y = Y.*mX;
%-----align AnaT by qxT, qyT----
    FAnaT = advect(hmemT(:,:,2),qxyT(:,:,1),qxyT(:,:,2));
    fileID5 = fopen('/scratch3/BMC/det/Jie.Feng/Jing/readoutmem/displace_vector_neg.dat','w')
  %  fileID6 = fopen('/scratch3/BMC/det/Jie.Feng/Jing/readoutmem/FA_field1.dat','w')
	FAna = transpose(FAnaT)
	fwrite(fileID5,qxy,precision)
%	fwrite(fileID6,FAna,precision)
	disp(qxy)

e = cputime-t

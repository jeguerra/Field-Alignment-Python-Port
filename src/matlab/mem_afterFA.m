%close all;
clear all;
t1 = clock
%Field alignment
%input:
%---two fields---
%A and B
%output:
%--C--> qx, qy : vector field
%--D align B to A


%=====input files====
%dir = ...    filename = ...
    fileID1 = fopen('/scratch3/BMC/det/Jie.Feng/Jing/readoutmem/GH_20mem_500hPa_0to360andS90toN90_120h.dat');
%====parameters====
    precision = 'float';
    nlon = 360; %dim of lon
    nlat = 181; %dim of lat
    nmem = 20;  %total number of ens members
    %----parameters for field alignment
    iter = 2^10;
    lscale = 1;   %change from 1~10,   1 is the best 

%====Initialize variables====
    qxyT_mem = zeros(nlat,nlon,2,nmem,'single');
    qxyT = zeros(nlat,nlon,2,'single');
    qxy = zeros(nlon,nlat,2,'single');
    hmemT_FA = zeros(nlat,nlon,'single');
    hmem_FA = zeros(nlon,nlat,'single');

%====read data====
    hmem0 = fread(fileID1,nlon*nlat*nmem,precision);   % read the file Yall(360*181*20)
    hmem = reshape(hmem0,[nlon,nlat,nmem]);
    %----transpose data before input to subroutine---- 
    hmemT = zeros(nlat,nlon,nmem,'single');
    for i = 1:1:nmem;
	hmemT(:,:,i) = transpose(hmem(:,:,i)); 
    end;

%=====Start calculation====
   imem=1;    %+++++This would be changed for parallel task+++
   fprintf('Start Calculation');
   fprintf('This is the member',imem)
%------derive the displacement vectors-------
for i = 1:1:nmem;
    fprintf('mem is:');
    [qxyT_mem(:,:,1,i),qxyT_mem(:,:,2,i)] = FA2DImNoHLIB(hmemT(:,:,imem),hmemT(:,:,i),128,.1,iter,2,lscale);
end;
%------calculate the mean of all DVs----
for i = 1:1:nmem;
    qxyT(:,:,1)=qxyT(:,:,1)+qxyT_mem(:,:,1,i);
    qxyT(:,:,2)=qxyT(:,:,2)+qxyT_mem(:,:,2,i);
end;
    qxyT(:,:,1)=qxyT(:,:,1)./nmem;
    qxyT(:,:,2)=qxyT(:,:,2)./nmem;
    qxy(:,:,1) = transpose(qxyT(:,:,1));    %output after transopose:
    qxy(:,:,2) = transpose(qxyT(:,:,2));    %output after transopose:
    fprintf('Derive the displacement vector successfully!')
%-----align the member i----
     hmemT_FA = advect(hmemT(:,:,imem),qxyT(:,:,1),qxyT(:,:,2));
     hmem_FA = transpose(hmemT_FA);    %output after transopose:
     fprintf('align member i successfully!')
%-----output data----
    filename1 = ['/scratch3/BMC/det/Jie.Feng/Jing/onetimecase/mem' int2str(imem) '_afterFA.dat'];
    filename2 = ['/scratch3/BMC/det/Jie.Feng/Jing/onetimecase/mem' int2str(imem) '_FAvector.dat'];
    fileID5 = fopen(filename1,'w');
    fileID6 = fopen(filename2,'w');
	fwrite(fileID5,hmem_FA,precision);
	fwrite(fileID6,qxy,precision);
t2 = clock
%e = etime(t2,t1)

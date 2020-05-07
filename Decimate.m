function [ decom,dim ] = Decimate( inpim )
%DECIMATE: Produce a 1-level Decomposition with a Laplacian Pyramid.
% Laplacian approximated by DoG. 

H = fspecial('gaussian',9,1);
om = filter2(H,inpim); %Low Pass
dim = inpim - om;  %High pass
decom = imresize(om,0.5,'bilinear'); %Decimated Low-pass
end


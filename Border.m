function newimage = Border(img, npix)
%% function newimage = Border(img, npix)
%% Copyright S. Chandu Ravela 1997-
%% This function adds a reflective border around an image. 
%% It is essential when filtering with large scales, for
%% otherwise the valid portion of the image becomes very small. 
%% This is used instead of zeroing out the border, doing circular filtering
%% or copying the last border pixel. All of those are biased. Reflection
%% leaves the filtered solution in a somewhat better state. 

sz = size(img);
newimage = zeros(sz+[2*npix 2*npix]);
newimage(npix+1:npix+sz(1),npix+1:npix+sz(2)) = img; %% Keep image intact in center
%% The remainder are reflection statements around the center.
newimage(npix+1:npix+sz(1),1:npix) = img(1:sz(1),npix:-1:1); 
newimage(npix+1:npix+sz(1),npix+sz(2)+1:2*npix+sz(2)) = img(1:sz(1),sz(2):-1:sz(2)-npix+1);
newimage(1:npix,npix+1:npix+sz(2)) = img(npix:-1:1,1:sz(2));
newimage(npix+sz(1)+1:2*npix+sz(1),npix+1:npix+sz(2)) = img(sz(1):-1:sz(1)-npix+1,1:sz(2));
newimage(1:npix,1:npix) = img(npix:-1:1,npix:-1:1);
newimage(npix+sz(1)+1:2*npix+sz(1),npix+sz(2)+1:2*npix+sz(2)) = img(sz(1):-1:sz(1) - npix +1,sz(2):-1:sz(2)-npix+1);
newimage(npix+sz(1)+1:2*npix+sz(1),npix:-1:1) = img(sz(1):-1:sz(1)-npix+1,1:npix);
newimage(npix:-1:1,sz(2)+npix+1:2*npix+sz(2)) = img(1:npix,sz(2):-1:sz(2)-npix+1);

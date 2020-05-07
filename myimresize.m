% interplate the field in(szin(1),szin(2)) to out(sz(1),sz(2))
% input:
%---in: the input array
%---sz: [rsz rsz]
% output:
%---out: the output array
function out = myimresize(in, sz)

szin = size(in);
xfrac = (szin(2)-1)/sz(2); yfrac = (szin(1)-1)/sz(1);
pitchx = 1+xfrac/2:xfrac:szin(2)-xfrac/2;
pitchy = 1+yfrac/2:yfrac:szin(1)-yfrac/2;
[x,y] = meshgrid(pitchx,pitchy);
[bdx,bdy] = meshgrid(1:szin(2),1:szin(1));
out = qinterp2(bdx,bdy,in,x,y,2);
%out = bicutest(bdx,bdy,in,x,y);

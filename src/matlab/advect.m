function XOut = advect(X0, qx, qy)
% Advect an image
%
% X0   - input image
% qx   - x vectors
% qy   - y vectors
% XOut - advected image
% Simple semi-lagrangian transport.

% If no vectors were given return input image
if(isempty(qx))
   XOut = X0;
   return
end
   qx(isnan(qx))=0;qy(isnan(qy))=0;
Xb = Border(X0,25);
[xx,yy] = meshgrid(1:size(X0,2),1:size(X0,1));

XOut = interp2(Xb,25+xx-qx,25+yy-qy, 'cubic',0);

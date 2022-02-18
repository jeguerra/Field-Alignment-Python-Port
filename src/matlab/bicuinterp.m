%%%bicubic is the user-defined function defined from: /usr/share/octave/2.9.9/m/general/bicubic.m

%% Copyright (C) 2005  Hoxide Ma
%
% This file is part of Octave.
%
% Octave is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% Octave is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, write to the Free
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.

% -*- texinfo -*-
% @deftypefn {Function File} {@var{zi}=} bicubic (@var{x}, @var{y}, @var{z}, @var{xi}, @var{yi})
%
% Return a matrix @var{zi} corresponding to the the bicubic
% interpolations at @var{xi} and @var{yi} of the data supplied
% as @var{x}, @var{y} and @var{z}. 
%
% For further information please see bicubic.pdf available at
% @url{http://wiki.woodpecker.org.cn/moin/Octave/Bicubic}
% @seealso{interp2}
% @end deftypefn

% Bicubic interpolation method.
% Author: Hoxide Ma <hoxide_dirac@yahoo.com.cn>

function F = bicubic (X, Y, Z, XI, YI, spline_alpha)

  if (nargin < 1 || nargin > 6)
    print_usage ();
  end

  if (nargin == 6 && prod (size (spline_alpha)) == 1)
    a = spline_alpha;
  else
    a = 0.5;
  end

  if (nargin <= 2)
    % bicubic (Z) or bicubic (Z, 2)
    if (nargin == 1) 
      n = 1;
    else
      n = Y;
    end
    Z = X;
    X = [];
    [rz, cz] = size (Z);
    s = linspace (1, cz, (cz-1)*pow2(n)+1);
    t = linspace (1, rz, (rz-1)*pow2(n)+1);
    
  elseif (nargin == 3)
      X
      Y
      
    if (~isvector (X) || ~isvector(Y))
      error ('XI and YI must be vector');
    end
    s = Y;
    t = Z;
    Z = X;
    [rz, cz] = size (Z);
  elseif (nargin == 5 || nargin == 6)
    [rz, cz] = size (Z) ; 
    
    if (isvector (X) && isvector (Y))
        size(X)
        size(Y)
        size(Z)
        rz
        cz
      if (rz ~= length (Y) || cz ~= length (X))
        error ('length of X and Y must match the size of Z');
      end
    elseif (size(X,1) == size(Y,1) && size(X,1) == size(Z,1))
      X = X(1,:);
      Y = Y(:,1);
    else
      error ('X, Y and Z must be martrices of same size');
    end
    
    % mark values outside the lookup table
    xfirst_ind = find (XI < X(1));
    xlast_ind  = find (XI > X(cz));    
    yfirst_ind = find (YI < Y(1));
    ylast_ind  = find (YI > Y(rz));
    % set value outside the table preliminary to min max index   
    XI(xfirst_ind) = X(1);
    XI(xlast_ind) = X(cz);
    YI(yfirst_ind) = Y(1);
    YI(ylast_ind) = Y(rz);


    X = reshape (X, 1, cz);
    ttmp = 1 + sign (X(cz))*eps;
    X(cz) = X(cz).*ttmp;
    if (X(cz) == 0) 
      X(cz) = eps;
    end; 
    size(XI) 
    XI = reshape (XI, 1, length(XI))
    [m, i] = sort ([X, XI]);
    o = cumsum (i <= cz);
    xidx = o(find (i > cz));
    
    Y = reshape (Y, rz, 1);
    ttmp = 1 + sign (Y(rz))*eps;
    Y(rz) = Y(rz).*ttmp;
    if (Y(rz) == 0) 
      Y(rz) = eps;
    end; 
    YI = reshape (YI, length (YI), 1);
    [m, i] = sort ([Y; YI]);
    o = cumsum (i <= rz);
    yidx = o([find( i> rz)]);
    
    % set s and t used follow codes
    s = xidx + ((XI - X(xidx))./(X(xidx+1) - X(xidx)));
    t = yidx + ((YI - Y(yidx))./(Y(yidx+1) - Y(yidx)));
  else
    print_usage ();
  end
  
  if (rz < 3 || cz < 3)
    error ('Z at least a 3 by 3 matrices');
  end

  inds = floor (s);
  d = find (s == cz);
  s = s - floor (s);
  inds(d) = cz-1;
  s(d) = 1.0;
  
  d = [];
  indt = floor (t);
  d = find (t == rz);
  t = t - floor (t);
  indt(d) = rz-1;
  t(d) = 1.0;
  d = [];

  p = zeros (size (Z) + 2);
  p(2:rz+1,2:cz+1) = Z;
  p(1,:)      =    (6*(1-a))*p(2,:)    -3*p(3,:)   + (6*a-2)*p(4,:);
  p(rz+2,:)   =    (6*(1-a))*p(rz+1,:) -3*p(rz,:)  + (6*a-2)*p(rz-1,:);
  p(:,1)      =    (6*(1-a))*p(:,2)    -3*p(:,3)   + (6*a-2)*p(:,4);
  p(:,cz+2)   =    (6*(1-a))*p(:,cz+1)  -3*p(:,cz) + (6*a-2)*p(:,cz-1);

  % calculte the C1(t) C2(t) C3(t) C4(t) and C1(s) C2(s) C3(s) C4(s)
  t2= t.*t;
  t3= t2.*t;

  ct0=     -a .* t3 +     (2 * a) .* t2 - a .* t ;      %# -a G0
  ct1 = (2-a) .* t3 +      (-3+a) .* t2          + 1 ;  %# F0 - a G1
  ct2 = (a-2) .* t3 + (-2 *a + 3) .* t2 + a .* t ;      %# F1 + a G0
  ct3 =     a .* t3 -           a .* t2;                %# a G1
  t = [];t2=[]; t3=[];

  s2= s.*s;
  s3= s2.*s;

  cs0=     -a .* s3 +     (2 * a) .* s2 - a .*s ;      %# -a G0
  cs1 = (2-a) .* s3 +    (-3 + a) .* s2         + 1 ;  %# F0 - a G1
  cs2 = (a-2) .* s3 + (-2 *a + 3) .* s2 + a .*s ;      %# F1 + a G0
  cs3 =     a .* s3 -           a .* s2;               %# a G1
  s=[] ; s2 = []; s3 = [];

  cs0 = cs0([1,1,1,1],:);
  cs1 = cs1([1,1,1,1],:);
  cs2 = cs2([1,1,1,1],:);
  cs3 = cs3([1,1,1,1],:);

  lent = length (ct0);
  lens = length (cs0);
  F = zeros (lent, lens);
  
  for i = 1:lent
    it = indt(i);
    int = [it, it+1, it+2, it+3];
    F(i,:) = [ct0(i),ct1(i),ct2(i),ct3(i)] * ...
        (p(int,inds) .* cs0 + p(int,inds+1) .* cs1 + ...
         p(int,inds+2) .* cs2 + p(int,inds+3) .* cs3);
  end

  % set points outside the table to NaN
  if (~ (isempty (xfirst_ind) && isempty (xlast_ind)))
    F(:, [xfirst_ind, xlast_ind]) = NaN;
  end
  if (~ (isempty (yfirst_ind) && isempty (ylast_ind)))
    F([yfirst_ind; ylast_ind], :) = NaN;
  end

  end

%!demo
%! A=[13,-1,12;5,4,3;1,6,2];
%! x=[0,1,4]+10; y=[-10,-9,-8];
%! xi=linspace(min(x),max(x),17);
%! yi=linspace(min(y),max(y),26);
%! mesh(xi,yi,bicubic(x,y,A,xi,yi));
%! [x,y] = meshgrid(x,y);
%! __gnuplot_raw__ ('set nohidden3d;\n')
%! hold on; plot3(x(:),y(:),A(:),'b*'); hold off;

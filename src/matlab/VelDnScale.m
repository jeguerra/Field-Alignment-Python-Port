function [ qqx,qqy ] = VelDnScale( qx,qy )
%A simple downscaling Vector-field resolution.

qqx = 2*imresize(qx,2);
qqy = 2*imresize(qy,2);
end


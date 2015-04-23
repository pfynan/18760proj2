function [ p ] = potential2( cx,cy,gx,gy,r )
%POTENTIAL2 Summary of this function goes here
%   Detailed explanation goes here

p = potential((abs(cx - gx)),r).*(potential((abs(cy - gy)),r));

end


function p = potential(d,r)
if (d>=0) && (d<=r/2)
    p = 1 - 2*(d^2)/(r^2);
elseif (d>=r/2) && (d<=r)
    p = 2*((d-r)^2)/(r^2);
elseif (d>r)
    p = 0;
end
end
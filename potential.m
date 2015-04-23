function p = potential(d,r)
p = ((d >= 0) & (d<=r/2)) .* (1 - 2*(d.^2)/(r^2));
p = p + ((d>r/2) & (d<=r)) .* (2*((d-r).^2)/(r^2));
end
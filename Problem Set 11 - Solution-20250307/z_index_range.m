function inds = z_index_range(zind, knum)
%z_index_range finds the indices for a given index of z level
firstind = (zind-1)*knum + 1;
lastind = (zind)*knum;
inds = firstind:lastind;
end
function ind = findMask(Xind,mask)

ind = find(all(Xind == mask*ones(1,size(Xind,2)),1));

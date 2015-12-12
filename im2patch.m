function patches = im2patch(I,patchSize,shiftSize)



if numel(size(I))<=2
patches = zeros(prod(patchSize),...
    floor((size(I,1)-patchSize(1)+1)/shiftSize(1))*floor((size(I,2)-patchSize(2)+1)/shiftSize(2)));
k = 0;

for j = 1:shiftSize(2):(size(I,2)-patchSize(2)+1)
    for i = 1:shiftSize(1):(size(I,1)-patchSize(1)+1)
        k = k+1;
        patches(:,k) = reshape(I(i:i+patchSize(1)-1,j:j+patchSize(2)-1),[prod(patchSize) 1]);
    end
end
else
    patches = zeros(prod(patchSize),...
        size(I,3)*floor((size(I,1)-patchSize(1)+1)/shiftSize(1))*floor((size(I,2)-patchSize(2)+1)/shiftSize(2)));
    k = 1-size(I,3);

    for j = 1:shiftSize(2):(size(I,2)-patchSize(2)+1)
        for i = 1:shiftSize(1):(size(I,1)-patchSize(1)+1)
            k = k+size(I,3);
            patches(:,k:k+size(I,3)-1) = reshape(I(i:i+patchSize(1)-1,j:j+patchSize(2)-1,:),[prod(patchSize) size(I,3)]);
        end
    end
end
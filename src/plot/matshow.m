function fh = matshow(img_mat)

fh = imagesc(img_mat);
axis equal
xlim([1, size(img_mat, 1)])
xlim([1, size(img_mat, 2)])
xlabel('u [px]')
ylabel('v [px]')

end


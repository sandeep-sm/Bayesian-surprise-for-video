figure(1);
for k = 1:1:nframes
    subplot(1,3,1);
    imshow(mat2gray(img(:,:,:,k)));
    title('Frame/scene');
    
    subplot(1,3,2);
    imshow(resimg(:,:,k));
    title('Saliency map');
    
    subplot(1,3,3);
    imshow(surprise(:,:,k));
    title('Surprise map');
    
    pause(0.2);
end

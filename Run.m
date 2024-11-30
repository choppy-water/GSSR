    original_image_path = '124084.jpg';
    marked_imageName = '124084_marked.png';
    mask1_path = '124084_mask1.png';
    mask2_path = '124084_mask2.png';
    
    Iorig = imread(original_image_path);
    Imarked = imread(marked_imageName);
    mask1 = imread(mask1_path);
    mask2 = imread(mask2_path);
  
    mask1(mask1 > 0) = 1;
    mask2(mask2 > 0) = 1;
    
    [m, n, ~] = size(Iorig);
    mask = zeros(m, n, 2);
    mask(:, :, 1) = mask1;
    mask(:, :, 2) = mask2;
    maskconstraints = mask;
    
    % Computing the Lp-based segmentation
    [~, Ibin] =  Lp_seg(Iorig, maskconstraints,0.1);
    res = res_deal(Ibin,  mask1);

    figure;
    imshow(res);
    % Cutting the segmentation
    Icut = LCcut(Iorig, res, 200);
    figure;
    imshow(mat2gray(Icut));
    
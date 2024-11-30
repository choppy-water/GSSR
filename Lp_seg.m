function [Isol, Ibin] =  Lp_seg(Iorig, maskconstraints, lambda)

% Computing the matrix of weights
disp('Computing the matrix of weights');
W = LCweights(Iorig);

% Building/solving the linear system of equations
disp('Building/solving the LC linear system');
[Lm,Lu,B, labels, cindex,index, xm, xsol] = L1linear_hard(W, maskconstraints); 

[Df,Db]=GMM(Iorig,maskconstraints);

rho=120 ;

d0=zeros(length(index),1);
b0=zeros(length(index),1);
d1=zeros(length(cindex),1);
b1=zeros(length(cindex),1);

w1_len = size(Lm * xm , 1);
w2_len = size(B' * xm , 1);

w1=spdiags(ones(w1_len,1),0,w1_len,w1_len);
w2=spdiags(ones(w2_len,1),0,w2_len,w2_len);

itr_num = 4;
for i=1:itr_num
    fprintf('Iteration %d out of %d...\n',i, itr_num);
    
    left_hand=(B'* w1^2 * B + Lu' * w2^2 * Lu) * rho + 2 * lambda * (Df + Db);
    right_hand = rho*(B'* w1^2 * Lm + Lu'* w2^2 * B') * xm ...
        -B'* w1 * (rho * d0 - b0) + Lu' * w2 * (rho * d1-b1) ...
        + 2 * lambda * Df * ones(size(Df,1),1);
    xunk = left_hand\right_hand;
    
    d0 = shrink(w1 * (Lm * xm - B * xunk) + b0 * 1.0 / rho,1.0 / rho);
    d1 = shrink(w2 * (Lu * xunk - B' * xm) + b1 * 1.0 / rho,1.0 / rho);
    
    b0 = b0 + w1 * (Lm * xm - B * xunk) - d0;
    b1 = b1 + w2 * (Lu * xunk - B' * xm) - d1;
    
    w1 = element_wise(Lm * xm - B * xunk , 0.8);
    w1 = spdiags(w1,0,w1_len,w1_len);
    w2 = element_wise(Lu * xunk - B'* xm , 0.8);
    w2 = spdiags(w2,0,w2_len,w2_len);
end

% Outlier pruning
disp('Segmenting the image');
xunk(xunk < labels(1)) = labels(1);
xunk(xunk > labels(2)) = labels(2);

% Assigning the obtained solution to the output vector
xsol(cindex) = xunk;
Isol = reshape(xsol, size(Iorig,1), size(Iorig,2));

% Generating the definitive segmentation
%op: 1 (by Otsu's thresholding)
Ibin = (Isol > graythresh(Isol));
%op: 2 (by trivial cutting)
% Ibin = (Isol > sum(labels)/2);

end


function d = shrink(v, lambda)
index = find(v < 0);
d = max(abs(v) - lambda, 0);
d(index) = -1 * d(index);
end

function W = element_wise(w , p)
w = abs(w);
w = w + 1e-5;
W= w.^(p-1);
end
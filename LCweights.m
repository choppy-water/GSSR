function W = LCweights(Iorig)

I = LCfeats(Iorig);

Epsilon = 0.0001;  
Beta = 600; 
[nRows, nCols, nChannels] = size(I);  

Img = zeros(nRows, nCols, nChannels); 
Dist_lr_aux = zeros(nRows-1, nCols-1, nChannels);
Dist_rl_aux = zeros(nRows-1, nCols-1, nChannels);

for k=1:nChannels  
   Img(:,:,k) = mat2gray(I(:,:,k));  
   Dist_lr_aux(:,:,k) = abs(Img(1:end-1,1:end-1,k) - Img(2:end,2:end,k));  
   Dist_rl_aux(:,:,k) = abs(Img(1:end-1,2:end,k) - Img(2:end,1:end-1,k));  
end

dist_lr = reshape(max(Dist_lr_aux, [], 3), (nRows-1)*(nCols-1), 1);
dist_rl = reshape(max(Dist_rl_aux, [], 3), (nRows-1)*(nCols-1), 1);
dist_ht = reshape(max(abs(diff(Img,1,1)), [], 3), (nRows-1)*nCols, 1);
dist_vt = reshape(max(abs(diff(Img,1,2)), [], 3), nRows*(nCols-1), 1);

dist = [dist_lr; dist_rl; dist_ht; dist_vt];  
dist = roundn(exp(-Beta*dist.*dist),-4) + Epsilon;

iaux = transpose(1:nRows*nCols);
   
ilr = iaux(1:end-nRows);
ilr(nRows:nRows:end) = [];  

jlr = ilr + nRows + 1;   

irl = ilr + 1;  
jrl = jlr - 1;  

iht = iaux;
iht(nRows:nRows:end) = []; 
jht = iht + 1;  

ivt = iaux(1:end-nRows);
jvt = iaux(nRows+1:end); 

i = [ilr; irl; iht; ivt];
j = [jlr; jrl; jht; jvt];

W = sparse([i; j], [j; i], [dist; dist], nRows*nCols, nRows*nCols, 2*numel(i));
end
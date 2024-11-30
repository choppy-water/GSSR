function Feats =  LCfeats(Iorig)

mm = ones(3)/(3)^2;  
rs = 5;

for i=1:size(Iorig,3)  
   u = conv2(Iorig(:,:,i), mm, 'same');  
   Iorig(:,:,i) = u;     
end

Feat1 = uint8(applycform(mat2gray(Iorig), makecform('srgb2lab')));
Feat2 = Iorig;

SA = rangefilt(Feat1, true(rs));   

Feats(:,:,1) = imadjust(uint8(SA(:,:,1)));
Feats(:,:,2) = imadjust(uint8(SA(:,:,2)));
Feats(:,:,3) = imadjust(uint8(SA(:,:,3)));

SB = rangefilt(Feat2, true(rs));
Feats(:,:,1+end) = imadjust(uint8(SB(:,:,1)));  
Feats(:,:,1+end) = imadjust(uint8(SB(:,:,2)));
Feats(:,:,1+end) = imadjust(uint8(SB(:,:,3)));
end
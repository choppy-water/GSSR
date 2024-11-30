function [Lm, Lu,B, labels, cindex,index,xm, xsol] = L1linear_hard(W, maskconstraints)

% Initializations
Eps = 1.0;
n = size(W,1);  
nRegions = size(maskconstraints,3);
Constrb = zeros(n,nRegions); 

labels = [0 1]; 
indc = 0;
for k=1:nRegions  
    mask = logical(maskconstraints(:,:,k));   
    indb = find(mask);   
    Constrb(indb,k) = 1;
end

DW = spdiags(sum(W,2), 0, n, n);    
L = Eps*(DW - W);  

b = labels(2)*Constrb(:,1) + labels(1)*Constrb(:,2);

index = 1:n;
index(~logical(Constrb(:,1) + Constrb(:,2))) = [];

cindex = 1:n;
cindex(index) = [];

Bt = L(cindex, index);  
B = -transpose(Bt);      
xm = b(index);         
Lu = L(cindex, cindex);
Lm = L(index, index);  

xsol = zeros(size(xm));  
xsol(index) = xm; 
end
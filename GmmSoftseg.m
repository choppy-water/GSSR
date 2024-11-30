classdef GmmSoftseg %< CoreBaseClass
   properties 
       FgGMM; 
       BgGMM;
       k;      % number of components
       iter;   % number of iteration
       lambda; % graph cut parameter, weight of pairwise potential (1.0-20.0)
       sigma;  % grapu cut parameter, control itensity difference (5.0-20.0)
   end
    methods
       function obj = GmmSoftseg()
           obj.lambda = 10.0;  
           obj.sigma = 40.0;
           obj.k = 3;
           obj.iter = 5;
           obj.FgGMM = gmdistribution();    
           
           obj.BgGMM = gmdistribution();    
       end
       function [Fprob,Bprob]=weights(obj,Iorig,mask)

           obj.lambda = 10.0;  
           obj.sigma = 40.0;
           obj.k = 5;
           obj.iter = 5;
           [H,W,C]=size(Iorig);

           Ireshape = double(reshape(Iorig,[H*W,C]));   
           Freshape = reshape(mask(:,:,1),[H*W,1]);

           Breshape = reshape(mask(:,:,2),[H*W,1]);
           BData = Ireshape(Breshape==1,:);

           Unknow=mask(:,:,1)+mask(:,:,2);
           UData=Ireshape(Unknow==0,:);
           options = statset('MaxIter',500);
           obj.FgGMM = fitgmdist(FData,obj.k, 'RegularizationValue',0.1,'Options',options);
           obj.BgGMM=  fitgmdist(BData,obj.k,'RegularizationValue',0.1,'Options',options);
           
           Fprob = obj.FgGMM.pdf(UData);
           Bprob = obj.BgGMM.pdf(UData);          
           
           i=size(Fprob);
           for m=1:i
               if(Fprob(m)>=Bprob(m))
                   Bprob(m)=0;
               else
                   Fprob(m)=0;
               end
           end
           
       end       
  end
end  


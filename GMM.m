function [wf,wb]=GMM(Iorig,mask)

unlabel_index=find(mask(:,:,1)+mask(:,:,2)==0);
m=size(unlabel_index,1);

gmm=GmmSoftseg();
[Fprob,Bprob]=gmm.weights(Iorig,mask);

 Beta=600;
 Epsilon = 0.0001;
 Fvalue = roundn(exp(-Beta*Fprob.*Fprob),-4) + Epsilon; 
 Bvalue=roundn(exp(-Beta*Bprob.*Bprob),-4) + Epsilon; 
 
 a=find(Fprob==0);
 Fvalue(a)=0;
  wf=spdiags(Fprob,0,m,m);
 
 b=find(Bprob==0);
 Bvalue(b)=0;

 wb=spdiags(Bprob,0,m,m);
end
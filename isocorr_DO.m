%distin  input isotopmer distribution matrix, one sample per column
%n: D number, m: O number
%im_D: D impurity, im_O: O impurity;
function [distout,mm,distout_pct]=isocorr_DO(distin,n,m,im_D,im_O)

ab_D=0.00001;  %natural 2D abundance
ab_O=0.00187; %natural 18O abundance

if nargin<4
   im_D=0.01;  %impurity 2D 
   im_O=0.01;  %impurity 18O
end

%D matrix
dim=(n+1)*(m+1);
M=zeros(dim,dim);
mm=zeros(dim,dim,4);
mmp=zeros(dim,dim,4);
for i=1:dim
   D_num_row=floor((i-1)/(m+1));
   O_num_row=mod((i-1),(m+1));
    for j=1:i
         D_num_col=floor((j-1)/(m+1));
         O_num_col=mod((j-1),(m+1));
         
         a1=n-D_num_col;
         b1=D_num_row-D_num_col;    
        
      
         a2=m-O_num_col;
         b2=O_num_row-O_num_col;      
      mm(i,j,:)=[a1,b1,a2,b2];
      if b1>=0 && b2>=0     
        M(i,j)=dbinom(a1,b1,ab_D)*dbinom(a2,b2,ab_O);
      end
    end
end

Mp=zeros(dim,dim);
for i=1:dim
   D_num_row=floor((i-1)/(m+1));
   O_num_row=mod((i-1),(m+1));
    for j=i:dim
         D_num_col=floor((j-1)/(m+1));
         O_num_col=mod((j-1),(m+1));
         
         a1=D_num_col;
         b1=D_num_col-D_num_row;    
        
      
         a2=O_num_col;
         b2=O_num_col-O_num_row;      
      mmp(i,j,:)=[a1,b1,a2,b2];
      if b1>=0 && b2>=0     
        Mp(i,j)=dbinom(a1,b1,im_D)*dbinom(a2,b2,im_O);
      end        

    end
end

N=M*Mp;
distout=distin;
 for i=1:size(distin,2)
     distout(:,i)=lsqnonneg(N,distin(:,i));%non-negative least square fitting
 end
 distout_pct=distout./(sum(distout,1)+1e-10);

% distout=max(distout,0); %nonzero
% distout=distout./sum(distout,2); 


function out=dbinom(a,b,r)
out=nchoosek(a,b)*(1-r)^(a-b)*r^b;

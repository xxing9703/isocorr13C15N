start_col=15;  %starting col number of data block.
impurity_C=0.01;
impurity_N=0.01;
%fname='20200310-M1-13C15N-Arg.csv';
fname='example1.csv';
%%
[filepath,name,ext] = fileparts(fname);
fname_S=fullfile(filepath,[name,'_cor',ext]);
A=readtable(fname);
compound=unique(A.compound,'stable');

dt_combo=[];
for i=1:length(compound)
 ids=find (strcmp(A.compound,compound(i)));
 A_sub=A(ids,:);
 
 [~,~,tp]=formula2mass(A_sub.formula{1});
 n=tp(1);m=tp(2);
 counts=[];
 lb=A_sub.isotopeLabel;
  for j=1:length(lb)
      str=lb{j};
      sub_str=split(str,'-');
      num=length(sub_str);
      if num==1
          counts(j,1)=0;
          counts(j,2)=0;          
      elseif num==3
          if strcmp(str(1),'C')
              counts(j,1)=str2num(sub_str{end});
              counts(j,2)=0; 
          elseif strcmp(str(1),'N')
              counts(j,1)=0;
              counts(j,2)=str2num(sub_str{end}); 
          else
              fprintf('something is wrong.......');              
          end
      elseif num==4
          counts(j,1)=str2num(sub_str{end-1});
          counts(j,2)=str2num(sub_str{end});          
      else
          fprintf('something is wrong');
      end
  end
  
 %aa=combvec(0:m,0:n)';
   v1=repmat(0:m,1,n+1);
   v2=reshape(repmat(0:n,m+1,1),1,(m+1)*(n+1));
   aa=[v1;v2]';
   
 
 aa=aa(:,[2,1]); %full list
 dt=A_sub{:,start_col:end};
 fulldt=[]; nr=[]; nj=[];
  for j=1:size(aa,1)
      tp=find(ismember(counts,aa(j,:),'rows'));
      if isempty(tp)
          fulldt=[fulldt;zeros(1,size(dt,2))];
      else
          nr=[nr,tp];
          nj=[nj,j];
          fulldt=[fulldt;dt(tp,:)];
      end
  end
 [~,~,dtout]=isocorr_CN(fulldt,n,m,impurity_C,impurity_N);
  [~,b]=sort(nr);
  dtout=dtout(nj(b),:);
  dt_combo=[dt_combo;dtout];
  
end
A_sub1=A(:,1:start_col-1);
A_sub2=A(:,start_col:end);
A_sub2{:,:}=dt_combo;
A_corrected=[A_sub1,A_sub2];

writetable(A_corrected,fname_S);





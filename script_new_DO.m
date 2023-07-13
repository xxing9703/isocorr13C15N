clear
fname='exp005_free_neg_omics.csv';
start_col=16;  %starting col number of data block.
impurity_D=0.01;
impurity_O=0.01;
autogrouping=0;
%%---------------------------------------------
A=readtable(fname,'readvariablename',true); 
A=A(1:length(find([A.medMz]>0)),:); %cut empty rows.
sample_name=A.Properties.VariableNames(start_col:end)';
grp_name=sample_name;
if autogrouping==1
     for i=1:length(sample_name)
         C=strsplit(sample_name{i},'_');
         if length(C)>1
             grp_name{i,1}=sample_name{i}(1:length(sample_name{i})-length(C{end})-1);
         else
             grp_name{i,1}=sample_name{i};
         end
     end 
end
grpHead=find(strcmp(A.isotopeLabel,'C12 PARENT'));
grpHead(end+1)=size(A,1)+1;
% read formula, extra D/O number   
      for i=1:length(grpHead)-1          
          ID(i)=i;          
          ids= grpHead(i):grpHead(i+1)-1;
          A.metaGroupId(ids)=ones(1,length(ids))*i;
        A_sub=A(ids,:); %A_sub: data sheet for the selected metabolite ID 
        try  %8/16/2020  added warning message for incorrect formulas
            [~,~,tp]=formula2mass(A_sub.formula{1});
            
        catch
            fprintf(['check row#',num2str(ids(1)+1),' for errors in the formula name: ',A_sub.formula{1}],'Error detected!');
            return
        end
        C_num=tp(3);  % D num
        N_num=tp(4);  % O num
     %short table-->full table (C+1)(N+1), insert rows with zero signals
        lb=A_sub.isotopeLabel; %string analysis to get number of C/N label
        counts=zeros(length(lb),2);
         for j=1:length(lb)
             str=lb{j};
             [Clb,Nlb,errmsg]=str2DO(str);
             if errmsg==0
               counts(j,1)=Clb;
               counts(j,2)=Nlb; 
             else
               fprintf('erros in isotopeLabel detected');
               return
             end
        end
     v1=repmat(0:N_num,1,C_num+1);
     v2=reshape(repmat(0:C_num,N_num+1,1),1,(N_num+1)*(C_num+1));
     cn=[v1;v2]';
     cn=cn(:,[2,1]); %full list,iterate all D/O combinations
     dt=A_sub{:,start_col:end};
     fulldt=[]; idx_r=[]; idx_j=[];
     for j=1:size(cn,1)
       tp=find(ismember(counts,cn(j,:),'rows'));
       if isempty(tp)
          fulldt=[fulldt;zeros(1,size(dt,2))];
       else
          idx_r=[idx_r,tp];
          idx_j=[idx_j,j];
          fulldt=[fulldt;dt(tp,:)];
       end
     end
       [~,b]=sort(idx_r);
       reduced_idx=idx_j(b);
        
        meta(i).ID=ID(i);      
        meta(i).name=A_sub.compound{1};
        meta(i).formula=A_sub.formula{1};
        meta(i).mz=A_sub.medMz(1);
        meta(i).rt=A_sub.medRt(1);
        meta(i).C_num=C_num;
        meta(i).N_num=N_num; 
        meta(i).original_abs=fulldt;
        meta(i).original_pct=fulldt./sum(fulldt,1);
        meta(i).tic=sum(fulldt,1);
        meta(i).reduced_idx=reduced_idx;        
      end
      
cat_abs=[];cat_pct=[];
for i=1:length(meta)  
 [corr_abs,~,corr_pct]=isocorr_DO(meta(i).original_abs,meta(i).C_num,meta(i).N_num,impurity_D,impurity_O);
 idx=meta(i).reduced_idx;  
 meta(i).corr_abs=corr_abs;
 meta(i).corr_pct=corr_pct;
 meta(i).corr_abs_short=corr_abs(idx,:); %shottable for output
 meta(i).corr_pct_short=corr_pct(idx,:); %shorttable for output
 meta(i).corr_tic=sum(meta(i).corr_abs,1);
  
 cat_abs=[cat_abs;corr_abs(idx,:)];  %concatenate for csv output
 cat_pct=[cat_pct;corr_pct(idx,:)];  
end

A_part1=A(:,1:start_col-1);
A_part2=A(:,start_col:end);

A_part2{:,:}=cat_abs;
A_corr_abs=[A_part1,A_part2];

A_part2{:,:}=cat_pct;
A_corr_pct=[A_part1,A_part2];

% make 3rd table as total ion/////////////
for i=1:length(meta)
  A_corr_total{i,1}=meta(i).ID;
  A_corr_total{i,2}=meta(i).name;
  A_corr_total{i,3}=meta(i).formula;
  for j=1:length(meta(i).tic)
  A_corr_total{i,3+j}=meta(i).tic(j);
  end
end
A_corr_total=cell2table(A_corr_total);
A_corr_total.Properties.VariableNames=[{'ID','Name','formula'},A.Properties.VariableNames(start_col:end)];
% end 3rd table  ////////////////////////

%save 
[filepath,name,~] = fileparts(fname);
% fname_S1=fullfile(filepath,[name,'_pct_cor','.csv']); %fname_S: filename for save
% fname_S2=fullfile(filepath,[name,'_abs_cor','.csv']);
% fname_S3=fullfile(filepath,[name,'_total_cor','.csv']);
% writetable(A_corr_pct,fname_S1);
% writetable(A_corr_abs,fname_S2);
% writetable(A_corr_total,fname_S3);

fname_all=fullfile(filepath,[name,'_cor','.xlsx']);
writetable(A,fname_all,'Sheet','original');
writetable(A_corr_pct,fname_all,'Sheet','cor_pct');
writetable(A_corr_abs,fname_all,'Sheet','cor_abs');
writetable(A_corr_total,fname_all,'Sheet','total');


fprintf('done!\n')

%% report
%plot2bar(meta,'corr_pct',1,grp_name,[6,6],'tmp');




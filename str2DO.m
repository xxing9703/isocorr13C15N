%This function parses string of isotopeLabel(from Maven output) and
%finds out C and N labels. Deliminator '-' is used to make the judgement.
%examples:  str='C13-label-1', [Cnum,Nnum]=[1,0]
%str='N15-label-1',[Cnum,Nnum]=[0,1]
%str='C13N15-label-2-1',[Cnum,Nnum]=[2,1]
%str='C12 PARENT', [Cnum,Num]=[0,0]

function [Dnum,Onum,errmsg]=str2DO(str)
      errmsg=0;
      sub_str=split(str,'-');
      num=length(sub_str);
      if num==1
          Dnum=0;
          Onum=0;          
      elseif num==3
          if strcmp(str(1),'D')              
              Dnum=str2num(sub_str{end});             
              Onum=0; 
          elseif strcmp(str(1),'O')
              Dnum=0;
              Onum=str2num(sub_str{end}); 
          else
              fprintf('it only works for D/O label'); 
              errmsg=1;
          end
      elseif num==4          
            Dnum=str2num(sub_str{end-1});
            Onum=str2num(sub_str{end});          
      else
          fprintf('something is wrong with the string');
          errmsg=2;          
      end
      if isempty(Dnum)||isempty(Onum)
          fprintf('labeling not a number');
          errmsg=3;          
      end

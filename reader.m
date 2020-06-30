%This is the reader input for zebrafish simulations

 prompt = 'Which simulation would you like to run? e.g. WT/nac/shd/pfe/rse/cho/seurat/leo/ablation/all : ';
 fish_type_str = input(prompt,'s');
%
 prompt = 'What will you name this simulation';
 name = input(prompt,'s');
%
%type_str={'shd'};
type_str={'shd_pfe', 'nac_pfe', 'nac_shd'};

 prompt = 'How many repeats would you like to do?  ';
 rep_str=input(prompt,'s');
%
 prompt = 'Would you like to record every dev stage, every half day or only final point (J+)? 1/2/3:';
 time_str=input(prompt,'s');
%
 str='Ok, I will generate %s repeat(s) of the fish type, %s' ;
 fprintf(str, rep_str, fish_type_str);
 
 for i=1:str2double(rep_str)
     if strcmp(fish_type_str,'all')==1
         for j=1:length(type_str)
             fish_type_str=type_str{j};
             name_dir=[fish_type_str,'_', num2str(i)];
             mkdir(name_dir);
             [domain_matrix,domain_matrix_ir,domain_matrix_X]=main(fish_type_str, str2double(time_str),name_dir);
             parsave(strcat(name_dir,'/',fish_type_str, '.mat'),domain_matrix,domain_matrix_ir,domain_matrix_X);
         end
     else
         name_dir=[name, fish_type_str, num2str(i)];
         mkdir(name_dir);
         [domain_matrix,domain_matrix_ir,domain_matrix_X]=main(fish_type_str, str2double(time_str),name_dir);
         parsave(strcat(name_dir,'/',fish_type_str, '.mat'),domain_matrix,domain_matrix_ir,domain_matrix_X);
         
     end
 end









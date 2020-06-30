function []= PROF_SAVE_NOW(time_str, mutant_type,domain_matrix, domain_matrix_ir, domain_matrix_X,rec,t,name_dir,SL)

stage_rec={'PB','PR','SP', 'SA','J', 'J+'};

if time_str==2
    %Every day save
parsave_2(strcat(name_dir,'/',mutant_type,'t=',num2str(round(t/(60*24))),'.mat'),t,SL,domain_matrix,domain_matrix_ir,domain_matrix_X);
else
    %Every dev stage save
parsave_2(strcat(name_dir,'/',mutant_type,'t=',num2str(round(t/(60*24))),'STAGE',stage_rec{rec},'.mat'),t,SL,domain_matrix,domain_matrix_ir,domain_matrix_X);
end
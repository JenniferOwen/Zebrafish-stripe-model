function [done_cho, count_id, count_il, count_i, count_xd, count_xb, domain_matrix_ir, domain_matrix_X, count_x]=check_timed_events(xan_on, mutant_type, count_id, count_il, count_i, count_xd, count_xb, domain_matrix_ir, done_cho, SL, size_rec,t, domain_matrix_X,listx)
%LIST REFERENCE
m=1;
xb=2;
xd=4;
il=5;
id=6;

[sizex, sizey]=size(domain_matrix_ir);
if strcmp(mutant_type,'cho')==1 && SL>=size_rec(2) && done_cho==0
    count_id=0;
    done_cho=1;
    while count_il<200
        R=randsample(sizex,1);
        C=randsample(sizey,1);
        if domain_matrix_ir(R,C)==0
            domain_matrix_ir(R,C)=il;
            count_il=count_il+1;
            
        end
    end
end

if t<1.5*60*24 && contains(mutant_type,'shd')==0 && ...
        strcmp(mutant_type,'ablation')==0 && strcmp(mutant_type,'ablate_iridophore')==0 && ...
        strcmp(mutant_type,'cho')==0 && strcmp(mutant_type,'rse')==0 && strcmp(mutant_type,'move_dwon_stripe')==0 ...
        && strcmp(mutant_type,'vertical_stripe')==0
     domain_matrix_ir(ceil(sizex/2)-1:ceil(sizex/2)+1,1:sizey)=6*ones(3,sizey);
    count_id=nnz(domain_matrix_ir==id);
    count_il=nnz(domain_matrix_ir==il);
    
elseif t<1.5*60*24 && strcmp(mutant_type,'move_down_stripe')==1
    domain_matrix_ir(ceil(sizex/4)-1:ceil(sizex/4)+1,1:sizey)=6*ones(3,sizey);
    count_id=nnz(domain_matrix_ir==id);
    count_il=nnz(domain_matrix_ir==il);
elseif t<1.5*60*24 && strcmp(mutant_type,'vertical_stripe')==1
    domain_matrix_ir(1:sizex,ceil(sizey/2)-1:ceil(sizey/2)+1)=6*ones(sizex,3);
    count_id=nnz(domain_matrix_ir==id);
    count_il=nnz(domain_matrix_ir==il);
elseif t>4.5*60*24 && t<5*60*24 && strcmp(mutant_type,'cho')==0 && strcmp(mutant_type,'ablation')==0 ...
         && xan_on==1 ...
        && strcmp(mutant_type,'move_down_stripe')==0  && strcmp(mutant_type,'vertical_stripe')==0 %delayed xanthophores in horizontal myoseptum
    domain_matrix_X(ceil(0.5*sizex)-1:ceil(0.5*sizex)+1,1:sizey)=xd*ones(3,sizey);
    count_xd=nnz(domain_matrix_X==xd);
    count_xb=nnz(domain_matrix_X==xb);
elseif t>4.5*60*24 && t<5*60*24 && strcmp(mutant_type,'move_down_stripe')==1 
    %delayed xanthophores in horizontal myoseptum
    domain_matrix_X(ceil(sizex/4)-1:ceil(sizex/4)+1,1:sizey)=xd*ones(3,sizey);
    count_xd=nnz(domain_matrix_X==xd);
    count_xb=nnz(domain_matrix_X==xb);
elseif t>4.5*60*24 && t<5*60*24 && strcmp(mutant_type,'vertical_stripe')==1 
    %delayed xanthophores in horizontal myoseptum
    %domain_matrix_X(3*ceil(sizex/4)-10:3*ceil(sizex/4)-8,1:sizey)=xd*ones(3,sizey);
    domain_matrix_X(1:sizex,ceil(sizey/2)-1:ceil(sizey/2)+1)=xd*ones(sizex,3);
    count_xd=nnz(domain_matrix_X==xd);
    count_xb=nnz(domain_matrix_X==xb);
elseif strcmp(mutant_type,'shd')==1
    domain_matrix_ir=zeros(sizex,sizey);
    count_id=0;
    count_il=0;
end

count_x=count_xd+count_xb;
count_i=count_id+count_il;
end

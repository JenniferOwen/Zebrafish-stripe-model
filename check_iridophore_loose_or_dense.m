function  [domain_matrix_ir,count_il,count_id] = check_iridophore_loose_or_dense(Q,D,m,count_il, count_id, xd,id, il, sizex,sizey,domain_matrix,domain_matrix_X,domain_matrix_ir,long_x,mutant_type)
%Function determines whether an iridophore in location
%iridophore_domain(R,C) transitions to opposing type or not.
 
%Inputs:
%[Q,D]: location of iridophore of interest
%m: reference for m-type
%count_il : number of loose iridophores on domain
%count_id : number of dense iridophores on domain
%xd: reference for xd
%id: reference for id
%il: reference for il
%[sizex,sizey]: size of iridophore domain
%domain_matrix: melanophore domain
%domain_matrix_X: xanthophore domain
%domain_matrix_I: iridophore domain
 
%Outputs:
%domain_matrix_ir: updated iridophore domain
%count_il: updated number of loose iridophores
%count_id: updated number of dense iridophores  
    
[neighbour_m]=calc_neighbours_xi_count_m(Q,D,m,2,sizex,sizey,domain_matrix,1);
[neighbour_LONG_x]=calc_neighbours(Q,D,xd,12,sizex,sizey,domain_matrix_X,1,0);
[neighbour_SHORT_x]=calc_neighbours(Q,D,xd,2,sizex,sizey,domain_matrix_X,1,1);


if strcmp(mutant_type,'leo_nac')==1 || strcmp(mutant_type, 'leo')==1  || strcmp(mutant_type, 'leo_sbr')==1
    if rand()<0.5
    neighbour_LONG_x=2; %bias loose to dense
    neighbour_SHORT_x=5;
    end
end

if strcmp(mutant_type,'luc')==1 
    if rand()<0.5
    neighbour_m=2; %bias loose to dense
    end
end


if domain_matrix_ir(Q,D)==id &&...
        ((neighbour_SHORT_x<9 && neighbour_LONG_x>14) || ... 
        (neighbour_m>4))
    %Transition from dense to loose
    domain_matrix_ir(Q,D)=il;
    count_il=count_il+1;
    count_id=count_id-1;
elseif domain_matrix_ir(Q,D)==il && ...
        ( neighbour_SHORT_x>4 || neighbour_LONG_x<16) ...
        && (neighbour_m<9 ) %12
    %Transition from loose to dense
    domain_matrix_ir(Q,D)=id;
    count_id=count_id+1;
    count_il=count_il-1; 
end


end

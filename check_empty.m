function [empty]=check_empty(R,C, domain_matrix,domain_matrix_X, domain_matrix_irr,m,sizex,sizey)
empty=1;

for i=-1:1
    for j=-1:1
        R_check=R+i;
        C_check=C+j;
        if (R_check>sizex || R_check<1)
            continue %periodic BCs and do not include the site that the cell is already occupying.
        end
        if C_check>sizey || C_check<0
            C_check=mod(C_check,sizey); %periodic BCs
        elseif C_check==0
            C_check=sizey; %periodic BCs
        end
        if m==0
            if (domain_matrix(ceil(R_check/2),ceil(C_check/2))==1 || domain_matrix_X(R_check,C_check)==4 || ...
                    domain_matrix_irr(R_check,C_check)==6)
                empty=0;
                break
            end
        else
            if (domain_matrix(R_check,C_check)==1 || nnz(domain_matrix_X(2*R_check-1:2*R_check,2*C_check-1:2*C_check)==4)>0 || ...
                    nnz(domain_matrix_irr(2*R_check-1:2*R_check,2*C_check-1:2*C_check)==6)>0)
                empty=0;
                break
            end
        end
        
    end
end
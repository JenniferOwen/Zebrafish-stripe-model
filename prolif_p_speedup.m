function [domain_matrix_1,addone]= prolif_p_speedup(R,C,domain_matrix_1, domain_matrix_X, domain_matrix, domain_matrix_ir, sizex, sizey)

%function takes cell on domain_matrix_1 in position R,C and executes
%proliferation. The new daughter cell is placed in one of four sites
%that is chosen at random if cell is a xanthoblast or iridophore. Or if a
%xanthophore it is dependent on the local distribution of cells.

event_chooser=rand();
addone=0;

%For xanthoblasts and iridophores, proliferation is random
if domain_matrix_1(R,C)==5 || domain_matrix_1(R,C)==6 || domain_matrix_1(R,C)==2
    if (event_chooser < 0.25)
        if (C ~= 1 ...
                && domain_matrix_1(R, C-1) == 0 )
            domain_matrix_1(R,C-1) = domain_matrix_1(R,C);
            addone=1;
        elseif C==1 && domain_matrix_1(R,sizey)==0 %periodic BCs
            domain_matrix_1(R,sizey) = domain_matrix_1(R,C);
            domain_matrix_1(R,C) = 0;
            addone=1;
        end
        %right proliferation
    elseif (event_chooser > 0.25 && event_chooser < 0.5)
        if (C ~= sizey ...
                && domain_matrix_1(R, C+1) == 0 )
            domain_matrix_1(R,C+1) = domain_matrix_1(R,C);
            addone=1;
        elseif C==sizey && domain_matrix_1(R,1)==0 %periodic BCs
            domain_matrix_1(R,1) = domain_matrix_1(R,C);
            domain_matrix_1(R,C) = 0;
            addone=1;
        end
        %up movement
    elseif (event_chooser > 0.5 && event_chooser < 0.75)
        if (R ~= 1 ...
                && domain_matrix_1(R-1, C) == 0 )
            domain_matrix_1(R-1,C) = domain_matrix_1(R,C);
            addone=1;
            
        end
        %down movement
    elseif (event_chooser > 0.75)
        if (R ~= sizex ...
                && domain_matrix_1(R+1, C) == 0 )
            domain_matrix_1(R+1,C) = domain_matrix_1(R,C);
            addone=1;
            
        end
    end
else %xanthophore proliferation
    
    %calculate the number of iridophores in short range
    [neighbour_id]=calc_neighbours(R,C,6,2,sizex,sizey,domain_matrix_ir,1,1);
    [neighbour_il]=calc_neighbours(R,C,5,2,sizex,sizey,domain_matrix_ir,1,1);
    [neighbour_xd]=calc_neighbours(R,C,4,2,sizex,sizey,domain_matrix_X,1,1);
    
    %calculate the number of melanophores
    [neighbour_m]=calc_neighbours_xi_count_m(R,C,1,2,sizex,sizey,domain_matrix,1);
    
    if (neighbour_id+neighbour_xd>2*(neighbour_il+neighbour_m))
        if (event_chooser < 0.25)
            if C~=1 && domain_matrix_1(R, C-1) == 0
                domain_matrix_1(R,C-1) = domain_matrix_1(R,C);
                addone=1;
            elseif C==1 && domain_matrix_1(R,sizey)==0 %periodic BCs
                domain_matrix_1(R,sizey) = domain_matrix_1(R,C);
                domain_matrix_1(R,C) = 0;
                addone=1;
            end
            %right movement
        elseif (event_chooser > 0.25 && event_chooser < 0.5)
            if C ~= sizey && domain_matrix_1(R, C+1) == 0
                domain_matrix_1(R,C+1) = domain_matrix_1(R,C);
                addone=1;
            elseif  C==sizey && domain_matrix_1(R,1)==0 %periodic BCs
                domain_matrix_1(R,1) = domain_matrix_1(R,C);
                domain_matrix_1(R,C) = 0;
                addone=1;
            end
            %up movement
        elseif (event_chooser > 0.5 && event_chooser < 0.75)
            if R ~= 1
                if domain_matrix_1(R-1, C) == 0
                    domain_matrix_1(R-1,C) = domain_matrix_1(R,C);
                    addone=1;
                    
                end
            end
            %down movement
        elseif (event_chooser > 0.75)
            if R ~= sizex
                if domain_matrix_1(R+1, C) == 0
                    domain_matrix_1(R+1,C) = domain_matrix_1(R,C);
                    addone=1;
                    
                end
            end
            
        end
    end
end
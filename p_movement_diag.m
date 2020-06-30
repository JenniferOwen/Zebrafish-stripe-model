function [domain_matrix]=p_movement_diag(domain_matrix, sizex,sizey , R, C, v)

%Normalised rates for each of the eight directions
A=(2/3)/(sum(v(1:4)));
B=(1/3)/(sum(v(5:8)));
v(1:4)=A*v(1:4);
v(5:8)=B*v(5:8);
v_prob=cumsum(v/sum(v));

%%Function describes random movement of a chosen particle in space into one
%%of eight directions with periodic BCs in the left and right and non-periodic up and down
%and volume exclusion from position (R,C)

movement_selector=rand();
%right movement
if (movement_selector < v_prob(1))
    if (C ~= sizey ...
            && domain_matrix(R, C+1) == 0)
        domain_matrix(R,C+1) = domain_matrix(R,C);
        domain_matrix(R,C) = 0;
    elseif C==sizey && domain_matrix(R,1)==0 %periodic BCs
        domain_matrix(R,1) = domain_matrix(R,C);
        domain_matrix(R,C) = 0;
    end
    %left movement
elseif (movement_selector > v_prob(1) && movement_selector < v_prob(2))
    if (C ~= 1 ...
            && domain_matrix(R, C-1) == 0)
        domain_matrix(R,C-1) = domain_matrix(R,C);
        domain_matrix(R,C) = 0;
    elseif C==1 && domain_matrix(R,sizey)==0 %periodic BCs
        domain_matrix(R,sizey) = domain_matrix(R,C);
        domain_matrix(R,C) = 0;
    end
    %up movement
elseif (movement_selector > v_prob(2) && movement_selector < v_prob(3))
    if (R ~= 1 ...
            && domain_matrix(R-1, C) == 0)
        domain_matrix(R-1,C) = domain_matrix(R,C);
        domain_matrix(R,C) = 0;
    end
    %down movement
elseif (movement_selector < v_prob(4) && movement_selector > v_prob(3))
    if (R ~= sizex ...
            && domain_matrix(R+1, C) == 0)
        domain_matrix(R+1,C) = domain_matrix(R,C);
        domain_matrix(R,C) = 0;
    end
    %Up and right movement
elseif (movement_selector < v_prob(5) && movement_selector > v_prob(4))
    if (C ~= sizey && R ~= 1  ...
            && domain_matrix(R-1, C+1) == 0)
        domain_matrix(R-1,C+1) = domain_matrix(R,C);
        domain_matrix(R,C) = 0;
    elseif C==sizey && R~=1 && domain_matrix(R,1)==0 %periodic BCs
        domain_matrix(R-1,1) = domain_matrix(R,C);
        domain_matrix(R,C) = 0;
    end
    %Up and left movement
elseif (movement_selector < v_prob(6) && movement_selector > v_prob(5))
    if (C ~= 1 && R~=1 ...
            && domain_matrix(R-1, C-1) == 0)
        domain_matrix(R-1,C-1) = domain_matrix(R,C);
        domain_matrix(R,C) = 0;
    elseif C==1 && R~=1 && domain_matrix(R,sizey)==0 %periodic BCs
        domain_matrix(R-1,sizey) = domain_matrix(R,C);
        domain_matrix(R,C) = 0;
    end
    %Down and right movement
elseif (movement_selector < v_prob(7) && movement_selector > v_prob(6))
    if (C ~= sizey && R~=sizex ...
            && domain_matrix(R+1, C+1) == 0)
        domain_matrix(R+1,C+1) = domain_matrix(R,C);
        domain_matrix(R,C) = 0;
    elseif C==sizey && R~=sizex && domain_matrix(R+1,1)==0 %periodic BCs
        domain_matrix(R+1,1) = domain_matrix(R,C);
        domain_matrix(R,C) = 0;
    end
    %Down and left movement
elseif movement_selector > v_prob(7)
    if (C ~= 1 && R~=sizex ...
            && domain_matrix(R+1, C-1) == 0)
        domain_matrix(R+1,C-1) = domain_matrix(R,C);
        domain_matrix(R,C) = 0;
    elseif C==1 && R~=sizex && domain_matrix(R+1,sizey)==0 %periodic BCs
        domain_matrix(R+1,sizey) = domain_matrix(R,C);
        domain_matrix(R,C) = 0;
    end
end
end





function [domain_matrix,domain_matrix_ir,domain_matrix_X]=main(mutant_type,time_str,name_dir)
close all

site_size_m=40;
day=24*60;

%LIST REFERENCE
m=1; %melanocytes
xb=2; %xanthoblasts
xd=4; %xanthophores
il=5; %loose iridophore
id=6; %dense iridophore

%INITIALISE NUMBERS OF CELLS ACCORDING TO MUTANT TYPE
[final_rec,N_xb, N_m, size_rec,xan_on,mel_on,irr_on]=determine_mutant(mutant_type);

%ALL PARAMETERS ARE INITIALISED HERE
[short_distance_m, short_distance_x, max_distance, long_x, growth_rateX,...
    growth_rateY, a_m, r_m, a_l, r_l, a_x, r_x, a_i, r_i, rand_move]=...
    fixed_parameters(mutant_type);

%INITIAL CONDITIONS FOR DOMAIN
[domain_matrix,domain_matrix_ir,domain_matrix_X,count_xb,count_m,count_id,count_xd,sizex_m,sizey_m,...
    sizex,sizey,~,count_i,count_il,t]...
    =Initial_conditions(mutant_type,N_xb,N_m,mel_on,irr_on);

%INITIALISE RECORDING TYPE
if time_str==1 %record every dev stage
    rec=1;
elseif time_str==2 %record every half day
    T=day;
    rec=1;
end

%Initialise whether timed events have occurred
done_cho=0;

SL=sizey_m*0.04+5.7;

while SL<=size_rec(final_rec) && t<60*day
    
    SL=sizey_m*0.04+5.7;
    
    %ALL TIMED EVENTS
    [done_cho, count_id, count_il, count_i, count_xd, count_xb, domain_matrix_ir, domain_matrix_X,count_x]=...
        check_timed_events(xan_on, mutant_type, count_id, count_il, count_i, count_xd, count_xb, ...
        domain_matrix_ir, done_cho, SL, size_rec,t, domain_matrix_X);
    
    
    %RECORDING THE DATA
    if time_str==1 %every dev stage
        if SL>=size_rec(rec)
            PROF_SAVE_NOW(1,mutant_type,domain_matrix, domain_matrix_ir, domain_matrix_X,rec,t,name_dir,SL);
            rec=rec+1; %store which stage we are in
        end
    elseif time_str==2 % every half day
        
        if t>=T
            PROF_SAVE_NOW(2,mutant_type,domain_matrix, domain_matrix_ir, domain_matrix_X,rec,t,name_dir,SL);
            T=T+60*24;
        end
        if SL>=size_rec(rec)
            rec=rec+1; %store which stage we are in
        end
    end
    
    
    %% DEFINE PARAMETERS
    [movement_m, diff_m, death_mel,...
        pull_xb, movement_il,  movement_xb,  movement_xd, movement_id,prolif_x, prolif_i...
        , ID_to_IL, IL_to_ID,diff_x]=...
        parameters(site_size_m, count_m, count_id, count_xd, count_xb,...
        count_i, count_il, count_x, mel_on,xan_on,irr_on,mutant_type);
    
    %% SUMMARISING ALL OF THE EVENT RATES
    diff=[diff_m, diff_x];
    movement_rate=[movement_m, movement_xb, movement_xd, movement_il, movement_id];
    prolif_rate=[prolif_x, prolif_i];
    transition_rates=[ID_to_IL, IL_to_ID];
    death_m=death_mel;
    pull=pull_xb;
    growth= [growth_rateX, growth_rateY];
    event_rates=[movement_rate, prolif_rate, diff, transition_rates, growth, death_m, pull];
    
    %GILLESPIE
    a0=sum(event_rates);
    a=cumsum(event_rates)/a0;
    tau  = (1/a0)*log(1/(rand()));
    t = t + tau;
    
    
    event_chooser=rand();
    
    if event_chooser<a(1) %movement_mel
        [R,C]=FIND(m,domain_matrix,sizex_m,sizey_m);
        [Attract,Repulse]=calc_AR_p_speedup_diag_m(R,C,xd,m,id,short_distance_m,sizex_m,sizey_m,...
            domain_matrix, domain_matrix_X, domain_matrix_ir);
        v=Attract*a_m+Repulse*r_m;
        if sum(v)==0
            v=rand_move;
        end
        [domain_matrix]=p_movement_diag(domain_matrix, sizex_m,sizey_m , R, C, v);
    elseif event_chooser<a(2) %movement_xb
        [R,C]=FIND(xb,domain_matrix_X,sizex,sizey);
        v=rand_move;
        [domain_matrix_X]=p_movement_diag(domain_matrix_X, sizex,sizey , R, C, v);
    elseif event_chooser<a(3) %movement_xd
        [R,C]=FIND(xd,domain_matrix_X,sizex,sizey);
        [Attract,Repulse]=calc_AR_p_speedup_diag(R,C,xd,m,id,short_distance_x,sizex,sizey,...
            domain_matrix, domain_matrix_X, domain_matrix_ir);
        v=Attract*a_x+Repulse*r_x;
        if sum(v)==0
            v=rand_move;
        end
        [domain_matrix_X]=p_movement_diag(domain_matrix_X, sizex,sizey , R, C, v);
    elseif event_chooser<a(4) %movement_il
        [R,C]=FIND(il,domain_matrix_ir,sizex,sizey);
        [Attract,Repulse]=calc_AR_p_speedup_diag(R,C,xd,m,il,short_distance_x,sizex,sizey,...
            domain_matrix, domain_matrix_X, domain_matrix_ir);
        v=Attract*a_l+Repulse*r_l;
        if sum(v)==0
            v=rand_move;
        end
        [domain_matrix_ir]=p_movement_diag(domain_matrix_ir, sizex,sizey , R, C, v);
    elseif event_chooser<a(5) %movement_id
        [R,C]=FIND(id,domain_matrix_ir,sizex,sizey);
%        if strcmp(mutant_type,'sbr')==0
            [Attract,Repulse]=calc_AR_p_speedup_diag(R,C,xd,m,il,max_distance,sizex,sizey,...
                domain_matrix, domain_matrix_X, domain_matrix_ir);
            v=Attract*a_i+Repulse*r_i;
            if sum(v)==0
                v=[1;1;1;1;1;1;1;1];
            end
%         elseif strcmp(mutant_type,'sbr')==1
%             v=[1;1;1;1;1;1;1;1];
%        end
        [domain_matrix_ir]=p_movement_diag(domain_matrix_ir, sizex,sizey , R, C, v);
    elseif event_chooser<a(6) %prolif x
        [R,C]=FIND_any(domain_matrix_X,sizex,sizey);
        [domain_matrix_X,addone]= prolif_p_speedup(R,C,domain_matrix_X, domain_matrix_X, domain_matrix, domain_matrix_ir, sizex, sizey);
        if addone>0 && domain_matrix_X(R,C)==xb
            count_xb=count_xb+1;
%            count_x=count_x+1;
        elseif addone>0 && domain_matrix_X(R,C)==xd
            count_xd=count_xd+1;
%            count_x=count_x+1;
        end
    elseif event_chooser<a(7) %prolif i
        [R,C]=FIND_any(domain_matrix_ir,sizex,sizey);
        [domain_matrix_ir,addone]= prolif_p_speedup(R,C,domain_matrix_ir, domain_matrix_X, domain_matrix, domain_matrix_ir, sizex, sizey);
        if addone>0 && domain_matrix_ir(R,C)==il
            count_il=count_il+1;
            count_i=count_i+1;
        elseif addone>0 && domain_matrix_ir(R,C)==id
            count_id=count_id+1;
            count_i=count_i+1;
        end
    elseif event_chooser<a(8) %mel differentiation
        R=randsample(sizex_m,1);
        C=randsample(sizey_m,1);
        
        %neighbour x_1, neighbour_il, neighbour_m are all inclusive of own space
        [neighbour_x_1]=calc_neighbours_m(R,C,xd,1,sizex_m,sizey_m,domain_matrix_X,1);
        [neighbour_m]=calc_neighbours(R,C,m,1,sizex_m,sizey_m,domain_matrix,4,1); %same cell type
        [neighbour_id]=calc_neighbours_m(R,C,id,1,sizex_m,sizey_m,domain_matrix_ir,1);
        
        [empty]=check_empty(R,C, domain_matrix,domain_matrix_X, domain_matrix_ir,1,sizex_m,sizey_m);
        
        %long m and long x are not inclusive of own space (add 4 m to
        %compare with x)
        [long_m]=calc_neighbours(R,C,m,6,sizex_m,sizey_m,domain_matrix,4,0); %same cell type
        [long_xd]=calc_neighbours_m(R,C,xd,6,sizex_m,sizey_m,domain_matrix_X,0);
        [long_id]=calc_neighbours_m(R,C,id,6,sizex_m,sizey_m,domain_matrix_ir,0);
        
        
        if contains(mutant_type,'leo')==0
            if((domain_matrix(R,C)==0) && (long_xd+long_id>2.5*long_m+3 && ... %3
                    neighbour_x_1<=neighbour_m && neighbour_id<=10) || ... %10
                    (rand()<0.01 && empty==1 ) )
                domain_matrix(R,C)=m;
                count_m=count_m+1;
            end %
        elseif (domain_matrix(R,C)==0) && (xan_on==0 || (long_xd>3*long_m && ... %3
                    neighbour_x_1<=neighbour_m)) %all leos
            domain_matrix(R,C)=m;
            count_m=count_m+1;
        end
    elseif event_chooser<a(9) %xanthoblast differentiation
        [R,C]=FIND(xb,domain_matrix_X,sizex,sizey);
        
        %calculate the number of iridophores in short range
        [neighbour_id]=calc_neighbours(R,C,id,2,sizex,sizey,domain_matrix_ir,1,1);
        [neighbour_il]=calc_neighbours(R,C,il,2,sizex,sizey,domain_matrix_ir,1,1);
        [neighbour_xd]=calc_neighbours(R,C,xd,2,sizex,sizey,domain_matrix_X,1,1);
        
        [empty]=check_empty(R,C, domain_matrix,domain_matrix_X, domain_matrix_X, 0,sizex,sizey);
        
        %calculate the number of melanophores in short range
        [neighbour_m]=calc_neighbours_xi_count_m(R,C,m,2,sizex,sizey,domain_matrix,1);
        
        if ((neighbour_id+neighbour_xd>(neighbour_il+neighbour_m))  || (rand<0.001 && empty==1))
            domain_matrix_X(R,C)=xd;
            count_xb=count_xb-1;
            count_xd=count_xd+1;
        elseif t>10*day && ((count_xd+count_id+count_il+count_m*4)/(sizex*sizey))<0.2 && (empty==1) && rand()<0.01
            domain_matrix_X(R,C)=xd;
            count_xb=count_xb-1;
            count_xd=count_xd+1;
        elseif t>30*day && ((count_xd+count_id+count_il+count_m*4)/(sizex*sizey))<0.4 && (empty==1) && rand()<0.01
            domain_matrix_X(R,C)=xd;
            count_xb=count_xb-1;
            count_xd=count_xd+1;
        end

    elseif event_chooser<a(10) %iridophore becomes loose (from presence of m)
        [R,C]=FIND(id,domain_matrix_ir,sizex,sizey);
        [domain_matrix_ir,count_il,count_id] = check_iridophore_loose_or_dense(R,C,m,count_il, count_id, xd,id, il, sizex,sizey,domain_matrix,domain_matrix_X,domain_matrix_ir,long_x,mutant_type);
    elseif event_chooser<a(11)
        [R,C]=FIND(il,domain_matrix_ir,sizex,sizey);
        [domain_matrix_ir,count_il,count_id] = check_iridophore_loose_or_dense(R,C,m,count_il, count_id, xd,id, il, sizex,sizey,domain_matrix,domain_matrix_X,domain_matrix_ir,long_x,mutant_type);
    elseif event_chooser<a(12) && ((sizex/sizey)<=1/2)
            % growth
        for j=1:2
            for i = 1:sizey
                aY_x=randsample(sizex,1);
                if rand()<0.5  %growth randomised
                    domain_matrix_ir((aY_x+1):(sizex+1),i) = domain_matrix_ir(aY_x:(sizex),i);
                    domain_matrix_ir(aY_x,i) = 0;
                else
                    domain_matrix_ir((aY_x+2):(sizex+1),i) = domain_matrix_ir(aY_x+1:(sizex),i);
                    domain_matrix_ir(aY_x+1,i) = 0;
                end
            end
            sizex=sizex+1;
        end
        sizex=sizex-2;
        
        for j=1:2
            for i = 1:sizey
                aY_x=randsample(sizex,1);
                if rand()<0.5
                    domain_matrix_X((aY_x+1):(sizex+1),i) = domain_matrix_X(aY_x:(sizex),i);
                    domain_matrix_X(aY_x,i) = 0;
                else
                    domain_matrix_X((aY_x+2):(sizex+1),i) = domain_matrix_X(aY_x+1:(sizex),i);
                    domain_matrix_X(aY_x+1,i) = 0;
                end
            end
            sizex=sizex+1;
        end
        
        
        for i = 1:sizey_m
            aY=randsample(sizex_m,1);
            if rand()>0.5  %centre the growth in the middle
                domain_matrix((aY+1):(sizex_m+1),i) = domain_matrix(aY:(sizex_m),i);
                domain_matrix(aY,i) = 0;
            else
                domain_matrix((aY+2):(sizex_m+1),i) = domain_matrix(aY+1:(sizex_m),i);
                domain_matrix(aY+1,i) = 0;
            end
        end
        sizex_m=sizex_m+1;
        %growth in DORSOLATERAL
    elseif event_chooser<a(13)
        for j=1:2
            for i = 1:sizex
                aX_x=randsample(sizey,1);
                if rand()>0.5  %centre the growth in the middle
                    domain_matrix_X(i,(aX_x+1):(sizey+1)) = domain_matrix_X(i,aX_x:(sizey));
                    domain_matrix_X(i,aX_x) = 0;
                else
                    domain_matrix_X(i,(aX_x+2):(sizey+1)) = domain_matrix_X(i,aX_x+1:(sizey));
                    domain_matrix_X(i,aX_x+1) = 0;
                end
            end
            sizey=sizey+1;
        end
        sizey=sizey-2;
        for j=1:2
            for i = 1:sizex
                aX_x=randsample(sizey,1);
                if rand()>0.5  %centre the growth in the middle
                    domain_matrix_ir(i,(aX_x+1):(sizey+1)) = domain_matrix_ir(i,aX_x:(sizey));
                    domain_matrix_ir(i,aX_x) = 0;
                    
                else
                    domain_matrix_ir(i,(aX_x+1):(sizey+1)) = domain_matrix_ir(i,aX_x:(sizey));
                    domain_matrix_ir(i,aX_x) = 0;
                end
            end
            sizey=sizey+1;
        end
        
        for i = 1:sizex_m
            aX = randsample(sizey_m,1);
            if rand()>0.5  %centre the growth in the middle
                domain_matrix(i,(aX+1):(sizey_m+1)) = domain_matrix(i,aX:(sizey_m));
                domain_matrix(i,aX) = 0;
                
            else
                domain_matrix(i,(aX+2):(sizey_m+1)) = domain_matrix(i,aX+1:(sizey_m));
                domain_matrix(i,aX+1) = 0;
                
            end
        end
        
        sizey_m=sizey_m+1;
    elseif event_chooser<a(14) % death of melanophores by xanthophores and melanophores
        [R,C]=FIND(m,domain_matrix,sizex_m,sizey_m);
        
        %neighbour x_1, neighbour_il, neighbour_m are all inclusive of own space
        [neighbour_x_1]=calc_neighbours_m(R,C,xd,1,sizex_m,sizey_m,domain_matrix_X,1);
        [neighbour_m]=calc_neighbours(R,C,m,1,sizex_m,sizey_m,domain_matrix,4,1); %same cell type
        [neighbour_il]=calc_neighbours_m(R,C,il,1,sizex_m,sizey_m,domain_matrix_ir,1);
        
        %long m and long x are not inclusive of own space. add 4 to m to
        %compare with x
        [long_m_6]=calc_neighbours(R,C,m,6,sizex_m,sizey_m,domain_matrix,4,0); %same cell type
        [long_xd]=calc_neighbours_m(R,C,xd,6,sizex_m,sizey_m,domain_matrix_X,0);
        
        if (strcmp(mutant_type,'leo_nac')==1 || strcmp(mutant_type, 'leo')==1 || ...
            strcmp(mutant_type, 'leo_shd')==1) || (strcmp(mutant_type, 'leo_pfe')==1 || strcmp(mutant_type, 'leo_sbr')==1) && ...
            (3*neighbour_m<neighbour_x_1)  || (rand<0.0001 && neighbour_il<10) %|| (rand<0.001 && neighbour_il<1))
            domain_matrix(R,C)=0;
            count_m=count_m-1;
        else
            if ((neighbour_m<neighbour_x_1) || (rand<0.001 && 3*long_xd<long_m_6 && neighbour_il<10))
            domain_matrix(R,C)=0;
            count_m=count_m-1;
            end
        end
    elseif event_chooser<a(15) %xanthoblast pulls melanophore
        [R,C]=FIND(xb,domain_matrix_X,sizex,sizey);
        if domain_matrix(ceil(R/2),ceil(C/2))==0
            [domain_matrix]=pull_xb_p(domain_matrix,R,C,sizex,sizey);
        end
    end
    
    
    
end

end






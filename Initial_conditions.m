function [domain_matrix,domain_matrix_ir,domain_matrix_X,count_xb,count_m,count_id,count_xd,sizex_m,...
    sizey_m,sizex_xi,sizey_xi,count_x,count_i,count_il,t...
    ]=Initial_conditions(mutant_type,N_xb,N_m,mel_on,irr_on)
%% INITIAL CONDITIONS

%WT, (seurat, sbr,leo): Iridophores along horizontal myoseptum (center), xanthoblasts randomly
%dispersed, melanophores dispersed randomly (Time = 0)

%ablation: A J+ stripe pattern with cells removed from a rectangular
%section. (Time = 30 days)

%ablate_iridophores: Iridophores along the horizontal myoseptum with a
%break in the middle. (Time = 0)

%rse: Same as WT with fewer iridophores (Time = 0)
%cho: No iridophores

%shd: No iridophores (Time = 0)
%pfe: No xanthophores (Time = 0)
%nac: No melanophores (Time = 0)


site_size_m=40; %lattice size on melanophore domain
site_size_xi=20; %lattice size on xanthophore/ iridophore domain.

if strcmp(mutant_type,'ablation')==1
    load('ablate.mat');

    %Ablate a section of cells from the domain
    domain_matrix(1:45,61:100)=zeros(45,40);
    domain_matrix_X(1:90,121:200)=zeros(90,80);
    domain_matrix_ir(1:90,121:200)=zeros(90,80);
    
    %Some interspersed iridophores and xanthoblasts remain
    for i=1:90
        for j=121:200
            if rand()<0.1
                domain_matrix_X(i,j)=2;
            end
            if rand()<0.1
                domain_matrix_ir(i,j)=0;
            else
                domain_matrix_ir(i,j)=5;
            end
        end
    end
    domain_matrix=domain_matrix(1:45,1:end);
    domain_matrix_ir=domain_matrix_ir(1:90,1:end);
    domain_matrix_X=domain_matrix_X(1:90,1:end);
    
    
    [sizex_xi,sizey_xi]=size(domain_matrix_X);
    [sizex_m,sizey_m]=size(domain_matrix);
    
    %Time is 30 days
    t=30*24*60;
    
elseif strcmp(mutant_type,'initially_stripey')==1
    sizex_m=floor(1000/site_size_m); %1000
    sizey_m=floor(2000/site_size_m);
    sizex_xi=floor(1000/site_size_xi); %1000
    sizey_xi=floor(2000/site_size_xi);
    
    %melanophore domain
    domain_matrix=zeros(sizex_m,sizey_m);
    %iridophore domain
    domain_matrix_ir=5*ones(sizex_xi,sizey_xi);
    %xanthophore domain
    domain_matrix_X=zeros(sizex_xi,sizey_xi);
    
    domain_matrix_ir(ceil(sizex_xi/2)-5:ceil(sizex_xi/2)+5,1:sizey_xi)=6*ones(11,sizey_xi);

    domain_matrix(2:9,1:sizey_m)=ones(8,sizey_m);
    domain_matrix(sizex_m-9:sizex_m-2,1:sizey_m)=ones(8,sizey_m);
  
    i=0;
     while i<N_xb %For all fish
        R=randsample(sizex_xi,1);
        C=randsample(sizey_xi,1);
        if domain_matrix_X(R,C)==0
            domain_matrix_X(R,C)=2;
            i=i+1;
        end
     end
    t=0;
    
    domain_matrix_X(ceil(sizex_xi/2)-5:ceil(sizex_xi/2)+5,1:sizey_xi)=4*ones(11,sizey_xi);
    domain_matrix_X(1:2,1:sizey_xi)=4*ones(2,sizey_xi);
    domain_matrix_X(sizex_xi-1:sizex_xi,1:sizey_xi)=4*ones(2,sizey_xi);
    domain_matrix_ir(1:2,1:sizey_xi)=6*ones(2,sizey_xi);
    domain_matrix_ir(sizex_xi-1:sizex_xi,1:sizey_xi)=6*ones(2,sizey_xi);
    
else %WT/ single-/ double- cell mutants/ rse
    
    
    %% Initialise domain
    %Initialise to be 2mm wide and 1mm height (at stage PB)
    sizex_m=floor(1000/site_size_m); %1000
    sizey_m=floor(2000/site_size_m);
    sizex_xi=floor(1000/site_size_xi); %1000
    sizey_xi=floor(2000/site_size_xi);
    
    if strcmp(mutant_type,'tall_domain')==1
        sizex_m=floor(3000/site_size_m); %1000
        sizey_m=floor(2000/site_size_m);
        sizex_xi=floor(3000/site_size_xi); %1000
        sizey_xi=floor(2000/site_size_xi);
    end
    
    if strcmp(mutant_type,'small_domain')==1
        sizex_m=floor(400/site_size_m); %1000
        sizey_m=floor(2000/site_size_m);
        sizex_xi=floor(400/site_size_xi); %1000
        sizey_xi=floor(2000/site_size_xi);
    end
    
    
    
    %iridophore domain
    domain_matrix_ir=zeros(sizex_xi,sizey_xi);
    %xanthophore domain
    domain_matrix_X=zeros(sizex_xi,sizey_xi);

    domain_matrix=zeros(sizex_m,sizey_m);
    if irr_on==1 %For all mutants/ WT fish with iridophores
        if strcmp(mutant_type,'rse')==1 %In the case of rse there are initially fewer iridophores
            r=1;
            while r<60
                R=randsample(3,1);
                C=randsample(sizey_xi,1);
                if domain_matrix_ir(R+ceil(sizex_xi/2)-2,C)==0
                    domain_matrix_ir(R+ceil(sizex_xi/2)-2,C)=6;
                    r=r+1;
                end
            end
        elseif strcmp(mutant_type,'rse')==0 && strcmp(mutant_type,'cho')==0 && strcmp(mutant_type,'move_up_stripe')==0 && strcmp(mutant_type,'vertical_stripe')==0            %Fish that are not rose or choker are initialised with a band of
            %iridophores along the horizontal myoseptum
            domain_matrix_ir(ceil(sizex_xi/2)-1:ceil(sizex_xi/2)+2,1:sizey_xi)=6*ones(4,sizey_xi);
        elseif strcmp(mutant_type,'move_down_stripe')==1
            %Fish that are not rose or choker are initialised with a band of
            %iridophores along the horizontal myoseptum
            domain_matrix_ir(ceil(sizex_xi/4)-1:ceil(sizex_xi/4)+1,1:sizey_xi)=6*ones(3,sizey_xi);
        elseif strcmp(mutant_type,'vertical_stripe')==1
            %Fish that are not rose or choker are initialised with a band of
            %iridophores along the horizontal myoseptum
            domain_matrix_ir(1:sizex_xi,ceil(sizey_xi/2)-1:ceil(sizey_xi/2)+1)=6*ones(sizex_xi,3);
        end
    end
    i=1;
    if mel_on==1 %For all mutants/ WT fish with melanophores
        while i<N_m
            %Add N_m melanophores unifomly at random to the domain
            R=randsample(sizex_m,1);
            C=randsample(sizey_m,1);
            if domain_matrix(R,C)==0
                domain_matrix(R,C)=1;
                i=i+1;
            end
        end
    end
    i=1;
    while i<N_xb %For all fish
        R=randsample(sizex_xi,1);
        C=randsample(sizey_xi,1);
        if domain_matrix_X(R,C)==0
            domain_matrix_X(R,C)=2;
            i=i+1;
        end
    end
    
    %Initialise time
    t=0;

        
    
end

%Initialise numbers of cells;
count_m=nnz(domain_matrix);
count_xb=nnz(domain_matrix_X==2);
count_xd=nnz(domain_matrix_X==4);
count_id=nnz(domain_matrix_ir==6);
count_x=count_xb; %no xanthophores
count_i=count_id;
count_il=0;

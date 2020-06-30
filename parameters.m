function [movement_m, diff_m, death_mel,...
    pull_xb, movement_il,  movement_xb,  movement_xd, movement_id,prolif_x, prolif_i...
    , ID_to_IL, IL_to_ID,diff_x]=...
    parameters(site_size_m, count_m, count_id, count_xd, count_xb,...
    count_i, count_il, count_x, mel_on,xan_on,irr_on,mutant_type)

%THESE PARAMETERS CHANGE WITH TIME
day=60*24;
site_size=20;

%Initialise parameters to be off
movement_m=0;
diff_m_i=0;
diff_m_x=0;
death_mel=0;
pull_xb=0;
diff_m_rand=0;
movement_il=0;
movement_id=0;
prolif_i=0;
ID_to_IL=0;
IL_to_ID=0;
movement_xb=0;
movement_xd=0;
prolif_x=0;
diff_x=0;

if mel_on==1
    %Attempted movement of melanocytes
    movement_m=count_m*110/(site_size_m*7*day);
    %Differentiation of m initiated by iridophores
    diff_m_i=count_id/(2*day);
    %Differentiation of m initiated by xanthophores
    diff_m_x=count_xd/(10*day);
    %Random differentiation of melanocytes
    diff_m_rand=300/day;
    if contains(mutant_type,'seurat')==1 
        diff_m_i=(1/10)*diff_m_i;
        diff_m_x=(1/10)*diff_m_x;
        diff_m_rand=(1/20)*diff_m_rand;
    end
    
    %Melanophore death attempted (affected in seurat)
    if contains(mutant_type,'seurat')==1
        death_mel=100*count_m/day;
    else
        death_mel=count_m/(3*day);
    end
    %Xanthoblast attaches to melanocyte
    pull_xb=count_xb/(100*day);
end

if irr_on==1
    %Iridophore movement
    movement_il=30*count_il/(site_size*60*24);
    movement_id=30*count_id/(site_size*60*24);
    %Iridophore proliferation (affected in rse)
    if strcmp(mutant_type,'rse')==1
        prolif_i=1.2*count_i/(5*day);
    else
        prolif_i=1.2*count_i/(day);
    end
    %Iridophore transition (affected in sbr)
    if contains(mutant_type,'sbr')==0 
        ID_to_IL=count_id/(day);
    else
        ID_to_IL=count_id/(10*day);
    end
    IL_to_ID=count_il/(day);
end

if xan_on==1
    %Movement of xanthophores
    movement_xb=count_xb*33/(site_size*7*day);
    movement_xd=count_xd*33/(site_size*7*day);
    %Proliferation of x
    prolif_x=(count_x)/(7*60*24);
    %Differentiation of xb is once every 5 days
    diff_x=count_xb/(5*day);
end


diff_m=sum([diff_m_i,diff_m_x,diff_m_rand]);

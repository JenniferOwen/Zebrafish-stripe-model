function [final_rec, N_xb,N_m,size_rec,xan_on,mel_on,irr_on]=determine_mutant(mutant_type)

%DEFAULT PARAMETERS

xan_on=1;
N_xb=2000;

mel_on=1;
N_m=50;

irr_on=1;

size_rec=[7.2,8.6,9.6,10.4,11,13];

%MUTANT CHANGES
if contains(mutant_type,'pfe')==1 
    xan_on=0;
    N_xb=0;
end
if contains(mutant_type,'shd')==1 
    irr_on=0;
end

if contains(mutant_type,'nac')==1 
    mel_on=0;
    N_m=0;
end

if contains(mutant_type,'rse')==1
    N_m=25;
end

if contains(mutant_type,'sbr')==1
    size_rec=[7.5,9,10.2,11.6];
end

if contains(mutant_type,'abla')==1
    size_rec=20;
end

final_rec=length(size_rec);

end



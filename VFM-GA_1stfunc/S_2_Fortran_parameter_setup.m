% ------------VFM GA hyperelastic V1.0---------------------------
% -------step2--the-fortran-codepack-governing-setup-------------
% ------------------------Zicheng Yan----------------------------
% -----Henann Research Group, Brown University-------------------
%%
clear all;
%this code deletes the existing, generate the new constants.f90 that
%controls the fortran code package under genetic algorithm optimization
%based on your choice of parameters

%% 1. Parameters that you should tune
run_mode = 2;
%0 for debug
%1 for special debug
%2 for formal run

alpha = 1.5; %2.15 
%read the paper, if lateral performance compared to axial is better, 
%reduce it, Or vice verse
%commonly take the range of 0.5 and 5, recommend between 1 and 3.5 for this
%model, best around 2-2.5 for first functionality
%is input data has defects (absense of lateral data),
%(input mode 3,input mode 4), code automatically adjust the alpha 

npop = 100;
%it determines how many independent populations you have, consider the
%parallelization, total 300 independent is the least requirement, if you
%parallelize by 100 batches, choose a number to be equal or larger than 3,
%Or you can simply have a number that is large, and shut down all fortran
%exe whenever you feel it's enough 

control_ps = 2; 
%2 for plane stress, which works the best, 1 for plane strain, which is not
%useful as most experiments are much closer to plane stress approximation,
%but plane strain can be verified through computational input data

lowb= [36447 33192 0.1 0.05 -0.0  0.05  2.0  0.01 0.05 2.00 2.00 5.00 0.001 1.0];
upb=  [154670 149789 0.3 6.0  -0.5  0.4  10.0  1.0  6.00 8.00 8.00 25.0 0.5   6.0];
% very important, input material parameter range that defines scope of
% optimiztion, if you are using the foam hyper model, no need to tune the
% 3rd to 14th parameter.

% for the first parameters G and K, estimate them by calculating the E and
% v in the small strain region (no code since data points are small strain
% could have a lot of corner cases), fill G and B estimated - 20% in lowb
% first 2 entries, fill estimated +20% in upb first 2 entries

% if the genetic algorithm code keep getting its best parameter set reaching only
% the maximum range --- upb, then you could increase the range for that parameter in upb


%These controls the use of Cmono and Cvol in 1st func
phys_c_e = 0;
phys_c_l = 0;
%0 for no use
%1 for Cmono = 1
%2 for Cvol = 1
%3 for Cmono = 1 & Cvol = 1 
% e for 1st gen check, l for later gen check
%only for 1st functionality unaxial experimental data, keep them the same
%if the code runs without problem, if you encounter that in the first 10
%generations, best(obj) is always 100000000, change both to 0, then
%gradually increase the parameter from 0 to 1, 0 to 2 and 0 to 3 while the
%code was able to produce actual value than 100000000
% actually this check is worse than stability check mentioned below, but
% still kept there



%% 2. Parameters that you can tune (you can also keep them (recommended))
ngen = 100; %number of generations per independent population
nt = 500; %number of parameter sets in a population, can be 300 can be 1000
xrate = 0.84; %cross over happen rate, 0.84 is great for 1st func
xcon = 0.9; %cross over rate for each single parameter in para set, 0.9 is great for 1st func
mrate = 0.7; %mutation rate happen rate
mcon = 0.2; %mutation rate for each single para in the set
mnum = 1; %1 is good as you only have 4 bits for each para in common setup

control_boltz = 1; 
%1 for boltzmann selection strongly recommened, 0 for roulette wheel
Tmax = 60;
Tmin = 0.5; 
%60 and 0.5 are the most common used boltzman temperature in boltzmann
%selection, you can change Tmax to 120 to allow more variable solution at
%first, it needs to be higher than the best objective function value
%computed in generation 0 ！！！！！！！！！！！！
cross_mode = 1; 
%1 for single point(recommended), 2 for shuffle, 3 for uniform 
mut_mode = 1;
%1 for simple mutation(recommend), 2 for rigged mutation 
elitism = 1;
%1 for elitism, 0 for deactivation, 1 strong recommend
uni_penalty = 0; 
%penalty for plane stress computing different out of plane deformation
%gradients in 1st functionality with uniaxial data, 1 for on, 0 for off

%(K1,K2)_c in paper
nstac = 0; %total number of compression stability check points
c_sta_k1=[-0.15,-0.25]%,-0.35,-0.45,-0.55]; %check points K1 value
c_sta_k2=[0.15,0.25]%,0.35,0.45,0.55]; %check points K2 value
% these 3 parameters defines the checks points at K3 = -1 for compression
% that checks the stability of material parameter set based on Dacorogna
% paper, check S_extra_ellipticity_Dacorogna.m for building a stability map
% You can add or reduce the check points in c_sta_k1 and c_sta_k2. If you
% see check points are unnecessary, reduce it, if you notice non monotonic
% behavior, add more check points, always play around the
% S_extra_ellipticity_Dacorogna.m before tuning it. 

%(K1,K2)_s in paper
s_sta = 0.51; %check point K2 value at K3 = 0 and K1 = 0 for simple shear
% similar to the 3 parameters above, this checks for shear stability, you
% can check  S_extra_ellipticity_Dacorogna.m for K3 = 0, simple shear stability


%% 3. parameters that you shouldn't tune 

di = 8; %precision of computation, 8 is the same for double in Abaqus, fortran default = 4
group = 4; %dependent on the sorting of external and interal force matrix, dont change it, mostly likely 4
sqnn = 52; %for indentation only, for computational verification, no longer useful, but still kept
bctype = 2; %no need to change for 1st functionality
intn = 4;%4 integration points used 
%1 for four face fictional Bc 
%2 mostly used, for simple compression/tension experiment
%3 for 3 face fictional Bc
%4 for indentaion verification computationally 

%% 4. parameters that you should tune if using this package on different constitutive model
nparam =length(lowb); %for https://github.com/HenannResearchGroup/ElastomericFoam
nbit = 4; %more number of bits represent larger number of potential parameters in the given range. 4 bits, 16 possibilities for a parameter, 16^14 enough for normal use

%% 5. parameters that not designed to tune
dim = 2; %dimension for data, is definitely 2D, as it's based on DIC with 2D approximation 
neln = 4; %usually 4 nodes per element, why not?

%the next parameters are read from S1 produced txt record file
fstruct = dir('*node.inp'); %dont put multiple node files in the folder, otherwise bug!
filename = fstruct.name;
nodename = fopen(filename,'r');
info_temp = textscan(nodename,'%d ,%f ,%f\n');
nnode = max(info_temp{1});
fstruct2 = dir('*element.inp');
filename = fstruct2.name;
elname = fopen(filename,'r');
info_temp = textscan(elname,'%d ,%d ,%d,%d ,%d\n');
nel = max(info_temp{1});
record_input_setup = fopen('record_input_data_setup.txt','r');
info_temp = textscan(record_input_setup,'%d %d %d\n');
data_mod = info_temp{1};
break_f = info_temp{2}+1;
nframes = info_temp{3}+info_temp{2};

if (data_mod == 1)
elseif (data_mod == 2)
    nframe = info_temp{2};
elseif (data_mod == 3) 
    alpha = 0.8;
elseif (data_mod == 4)
    nframe = info_temp{2};
    alpha = 0.8;
end

%%
location = [pwd '\GA_fort'];
loc = [location '\constants.f90'];
if isfile(loc)
    delete(loc)
end
filename = loc;
fid = fopen(filename,'w');
fprintf(fid,'module constants\n');
fprintf(fid,'implicit none\n');
fprintf(fid,'integer dim,neln,nparam,intn,nnode,nel,nframe,di,group,npop,lam_length\n');
fprintf(fid,'integer sqnn, break\n');
fprintf(fid,'real(%d) alpha\n',di);
fprintf(fid,'real(%d) lowb(%d),upb(%d),s_sta,c_sta_k1(%d),c_sta_k2(%d)\n',di,nparam,nparam,nstac,nstac);
fprintf(fid,'real(%d) xrate,mrate,xcon,mcon,keep_prob,Tmax,Tmin\n',di);
fprintf(fid,'integer nbit,ngen,nt,mnum,nstac\n');
fprintf(fid,'integer mut_mode,cross_mode,uni_penalty\n');
fprintf(fid,'integer control_boltz,elitism\n');
fprintf(fid,'integer phys_c_e,phys_c_l,control_ps,bctype,run_mode\n');

fprintf(fid,'parameter(dim=%d,neln=%d,nparam=%d,intn=%d,nnode=%d,nel=%d)\n',dim,neln,nparam,intn,nnode,nel);
fprintf(fid,'parameter(nframe=%d,di=%d,group=%d,npop=%d,sqnn=%d)\n',nframes,di,group,npop,sqnn);
fprintf(fid,'parameter(alpha=%3.3fd0,break=%d,nbit=%d,ngen=%d,nt=%d)\n',alpha,break_f,nbit,ngen,nt);
fprintf(fid,'parameter(xrate=%3.3fd0,mrate=%3.3fd0,xcon=%3.3fd0,mcon=%3.3fd0,mnum=%d)\n',xrate,mrate,xcon,mcon,mnum);
fprintf(fid,'parameter(nstac = %d, s_sta = %3.3fd0)\n',nstac,s_sta);
if nstac > 0
sta_format1 = 'parameter(c_sta_k1=(/%3.3fd0';
sta_format2 = 'parameter(c_sta_k2=(/%3.3fd0';
for i = 1:nstac-1
    sta_format1 = [sta_format1 ',%3.3fd0'];
    sta_format2 = [sta_format2 ',%3.3fd0'];
end
sta_format1 = [sta_format1 '/))\n'];
sta_format2 = [sta_format2 '/))\n'];
fprintf(fid,sta_format1,c_sta_k1);
fprintf(fid,sta_format2,c_sta_k2);
end
fprintf(fid,'parameter(mut_mode = %d ,cross_mode = %d,uni_penalty = %d)\n',mut_mode,cross_mode,uni_penalty);
fprintf(fid,'parameter(control_boltz = %d ,phys_c_e = %d, phys_c_l = %d)\n',control_boltz,phys_c_e,phys_c_l);
fprintf(fid,'parameter(Tmax = %3.3fd0,Tmin = %3.3fd0)\n',Tmax,Tmin);
fprintf(fid,'parameter(elitism = %d,control_ps = %d)\n',elitism,control_ps);
fprintf(fid,'parameter(bctype = %d)\n',bctype);
fprintf(fid,'parameter(run_mode = %d)\n',run_mode);

paralines = round(nparam/8); %print 8 paras per line
paraformat1 = 'parameter(lowb=(/%3.3fd0';
paraformat2 = 'parameter(upb=(/%3.3fd0'; 
line_count = 1;
flagg = 0;
for i = 1:nparam-1
    if (flagg == 0)
        paraformat1 = [paraformat1 ',%3.3fd0'];
        paraformat2 = [paraformat2 ',%3.3fd0'];
    elseif (flagg == 1)
        paraformat1 = [paraformat1 '%3.3fd0'];
        paraformat2 = [paraformat2 '%3.3fd0'];
        flagg = 0;
    end
    
    if i/(line_count*8) == 1 
        line_count = line_count+1;
        paraformat1 = [paraformat1 ', &\n &    '];
        paraformat2 = [paraformat2 ', &\n &    '];
        flagg = 1;
    end
end
paraformat1 = [paraformat1 '/))\n'];
paraformat2 = [paraformat2 '/))\n'];

fprintf(fid,paraformat1,lowb);
fprintf(fid,paraformat2,upb);
fprintf(fid,'end module constants');
fclose('all')


% ---------VFM GA hyperelastic V1.0-----------------------------
% ------------------2nd functionality---------------------------
% -------step2--the-fortran-codepack-governing-setup------------
% ---------------------Zicheng Yan------------------------------
% -----Henann Research Group, Brown University------------------


%%
clear all;
%this code deletes the existing, generate the new constants.f90 that
%controls the fortran code package under genetic algorithm optimization
%based on your choice of parameters, Go through at least part 1. part 2.

%% 1. Parameters that you should tune
run_mode = 2;
%0 for debug
%2 for formal run

num_mode = 2;
%numerical integration mod
% 1 for 4 integ points bisection+Newton Raphson (faster)
% 2 for 1 integ with reduced quadrature Bisection (slower, slightly more accurate)

nconfig = 2;%make sure all the required files are placed at /input_data folder
%if you are fitting 2 separate DIC experimental data together, set it to 2,
%otherwise 1. 
%When setting it to 1, prepare comp_node.inp, comp_element.inp
%and w00.rpt - w??.rpt according to instructions in
%S0_input_data_structuring.m

%When setting it to 2, prepare ten_node.inp, ten_element.inp and t00.rpt -
%t??.rpt according to instruction file.

alpha = 1.4;
%read the paper, if fitted lateral performance compared to axial is better, 
%reduce it, Or vice verse
%For second functionality, 1.3-2.0 is the recommended range

tweight = 0.7;
%This is the weight of the objective function contributed by the second configuration,
%e.g. fitting field
%data from simple compression and tension, since compression is harder to
%fit with complex response, set tweight below 1 to reduce the tension
%weight on fitting. 

ncomp = 17;%number of w??.rpt you have
nten = 13;%number of t??.rpt you have, put 0 if nconfig = 1

npop = 100;
%it determines how many independent populations you have, consider the
%parallelization, total 300 independent is the least requirement, if you
%parallelize by 100 batches, choose a number to be equal or larger than 3,
%Or you can simply have a number that is large, and shut down all fortran
%exe whenever you feel it's enough 


lowb= [82000 153800 0.1 0.05 -0.0  0.05  2.0  0.01 0.05 2.00 2.00 5.00 0.001 1.0];
upb=  [122000 233800 0.3 6.0  -0.5  0.4  10.0  1.0  6.00 8.00 8.00 25.0 0.5   6.0];
% very important, input material parameter range that defines scope of
% optimiztion, if you are using the foam hyper model, no need to tune the
% 3rd to 14th parameter.

% for the first parameters G and K, estimate them by calculating the E and
% v in the small strain region (no code since data points are small strain
% could have a lot of corner cases), fill G and B estimated - 20% in lowb
% first 2 entries, fill estimated +20% in upb first 2 entries

% if the genetic algorithm code keep getting its best parameter set reaching only
% the maximum range --- upb, then you could increase the range for that parameter in upb


phys_c_e = 1;
phys_c_l = 3; 
%0 for no use
%1 for Cmono = 1
%2 for Cvol = 1
%3 for Cmono = 1 & Cvol = 1 
%!!!!! This 2 variables need to be set to non zero ONLY when you are fitting to
%field data that are close to homogeneous (near perfect simple
%compression/tension experimental data)
%!!!!! If you are fitting to inhomogeneous data, set them to 0

%important because there could be 'bug' like conditions related to these
%two parameters, these are for 'rough' check for 2nd functionality

%0 for no rough check 
%1 for monotonic check
%2 for monotonic check
%3 for both checks 
%if you encounter that in the first 10
%generations, best(obj) is always 100000000, change both to 0, then
%gradually increase the parameter from 0 to 1, 0 to 2 and 0 to 3 while the
%code was able to produce actual value than 100000000
% actually this check is worse than stability check mentioned below, but
% they still work



%% 2. Parameters that you can tune (you can keep most of them (recommended))
ngen = 80; %number of generations per independent population
nt = 300; %number of parameter sets in a population, can be 300 can be 1000
xrate = 0.84; %cross over happen rate, 0.84 is great for 1st func
xcon = 0.9; %cross over rate for each single parameter in para set, 0.9 is great for 1st func
mrate = 0.03; %mutation rate happen rate
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
%1 for Elitism, 0 for deactivation, 1 strong recommend
uni_penalty = 0; 
%penalty for plane stress computing different out of plane deformation
%gradients in 1st functionality with uniaxial data, 1 for on, 0 for off
%not useful here


%(K1,K2)_c in paper
nstac = 5; %total number of compression stability check points
c_sta_k1=[-0.15,-0.25,-0.35,-0.45,-0.55]; %check points K1 value
c_sta_k2=[0.15,0.25,0.35,0.45,0.55]; %check points K2 value
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
bctype = 2; %no need to change for 
if num_mode == 1
    intn = 4;
    rflag = 0;
elseif num_mode == 2
    intn = 1;
    rflag = 1;
end
    

%1 for four face fictional Bc 
%2 mostly used, for Instron experiment, like simple compression/tension, or
%read the inhomogeneous case in my paper
%3 for 3 face fictional Bc
%4 for indentaion verification computationally 

%% 4. parameters that you should tune if using this package on different constitutive model
nparam =length(lowb); %for https://github.com/HenannResearchGroup/ElastomericFoam
nbit = 4; %more number of bits represent larger number of potential parameters in the given range. 4 bits, 16 possibilities for a parameter, 16^14 enough

%% 5. parameters that not designed to tune
dim = 2; %dimension for data, is definitely 2D, as it's based on DIC with 2D approximation 
neln = 4; %usually 4 nodes per element, why not?

cd input_data
fstruct = dir('comp_node.inp'); 
filename = fstruct.name;
nodename = fopen(filename,'r');
info_temp = textscan(nodename,'%d ,%f ,%f\n');
nnode = max(info_temp{1});
fstruct2 = dir('comp_element.inp');
filename = fstruct2.name;
elname = fopen(filename,'r');
info_temp = textscan(elname,'%d ,%d ,%d,%d ,%d\n');
nel = max(info_temp{1});

if nconfig == 2
    fstruct = dir('ten_node.inp'); 
    filename = fstruct.name;
    nodename = fopen(filename,'r');
    info_temp = textscan(nodename,'%d ,%f ,%f\n');
    nnode2 = max(info_temp{1});
    fstruct2 = dir('ten_element.inp');
    filename = fstruct2.name;
    elname = fopen(filename,'r');
    info_temp = textscan(elname,'%d ,%d ,%d,%d ,%d\n');
    nel2 = max(info_temp{1});
end

nframe = nten+ncomp;
nbreak = ncomp+1;
cd ..

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
fprintf(fid,'real(%d) tweight,alpha\n',di);
fprintf(fid,'real(%d) lowb(%d),upb(%d),s_sta,c_sta_k1(%d),c_sta_k2(%d)\n',di,nparam,nparam,nstac,nstac);
fprintf(fid,'real(%d) xrate,mrate,xcon,mcon,keep_prob,Tmax,Tmin\n',di);
fprintf(fid,'integer nbit,ngen,nt,mnum,nstac\n');
fprintf(fid,'integer mut_mode,cross_mode,uni_penalty\n');
fprintf(fid,'integer control_boltz,elitism\n');
fprintf(fid,'integer phys_c_e,phys_c_l,bctype,run_mode\n');
fprintf(fid,'integer ncomp,nten\n');
fprintf(fid,'integer rflag,nnode2,nel2\n');

fprintf(fid,'parameter(dim=%d,neln=%d,nparam=%d,intn=%d,nnode=%d,nel=%d)\n',dim,neln,nparam,intn,nnode,nel);
fprintf(fid,'parameter(nnode2=%d,nel2=%d,ncomp=%d,nten=%d,break=%d)\n',nnode2,nel2,ncomp,nten,nbreak);
fprintf(fid,'parameter(nframe=%d,di=%d,group=%d,npop=%d,sqnn=%d)\n',nframe,di,group,npop,sqnn);
fprintf(fid,'parameter(alpha=%3.3fd0,nbit=%d,ngen=%d,nt=%d)\n',alpha,nbit,ngen,nt);
fprintf(fid,'parameter(xrate=%3.3fd0,mrate=%3.3fd0,xcon=%3.3fd0,mcon=%3.3fd0,mnum=%d)\n',xrate,mrate,xcon,mcon,mnum);
fprintf(fid,'parameter(nstac = %d, s_sta = %3.3fd0,tweight = %3.3fd0)\n',nstac,s_sta,tweight);
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
fprintf(fid,'parameter(elitism = %d)\n',elitism);
fprintf(fid,'parameter(bctype = %d)\n',bctype);
fprintf(fid,'parameter(run_mode = %d)\n',run_mode);
fprintf(fid,'parameter(rflag = %d)\n',rflag);

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


% ---------VFM GA hyperelastic V1.0------------------------------
% ------------------2nd functionality----------------------------
% --------------step3--the-driver--------------------------------
% -----------------------Zicheng Yan-----------------------------
% -----Henann Research Group, Brown University-------------------


%this code MUST be run with 'Run Section' based on user's need

%see 1 if you wanna run the package on your own computer (windows based)
%(requires Intel oneAPI fortran compiler and visual studio) 
% https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.oj1rhe

%see 2 if you wanna run the package on cluster 
%(requires the cluster to have intel fortran compiler installed 
%(oneAPI or parallel studio XE)


%see 3 after the run and collect the 'best' material parameters
%plot it against the experimental data input
%% 1.1-----------run locally_starts_here----------------
%IMPORTANT NOTICE, run the code in where the code lies, dont change the
%working folder of MATLAB

%    COMPILE (manual)
%1. in the windows start menu, search for 'Intel oneAPI command prompt'
%2. open the command prompt manually
%3. type in: cd [the address of folder GA_fort] 
%4. copy the following and type in: 
%   ifort constants.f90 linalg.f90 fem_basics.f90 fem_const.f90 fem_compute_el.f90 stability.f90 fem_compute_global.f90 fem_preprocess.f90 print_write.f90 genetic_algo.f90 main.f90 -o run_ga.exe  

%OR

% COMPILE (AUTO THROUGH VBS)
%1. find the Intel one API command prompt location and name of that compiler
%and paste to one_API_address

oneAPI_address = 'C:\ProgramData\Microsoft\Windows\Start Menu\Programs\Intel oneAPI 2021\Intel oneAPI command prompt for Intel 64 for Visual Studio 2019';
%2. click Run Section
oneAPI = insertAfter(oneAPI_address,'\','\');
fortran_address = [pwd '\GA_fort'];
fortran = insertAfter(fortran_address,'\','\');
filename = 'compile.vbs';
fid = fopen(filename,'w');
str1 = ['oShell.Run("""',oneAPI,'""")\n'];
str2 = ['oShell.SendKeys "cd ',fortran,'{ENTER}"\n'];

fprintf(fid,'Set oShell = CreateObject("WScript.Shell")\n');
fprintf(fid,str1);
fprintf(fid,'WScript.Sleep 6000\n');
fprintf(fid,str2);
fprintf(fid,'WScript.Sleep 1000\n');
fprintf(fid,'oShell.SendKeys "ifort constants.f90 linalg.f90 fem_basics.f90 fem_const.f90 fem_compute_el.f90 stability.f90 fem_compute_global.f90 fem_preprocess.f90 print_write.f90 genetic_algo.f90 main.f90 -o run_ga.exe  {ENTER}"');
fclose('all');
system('compile.vbs')
%Dont interrupt(by hitting anything on computer) during VBS running

%% 1.2-----------Clear old data Create folders, copy files and run----------------
% change to the number of external parallelzation as you want, and hit run
% section, dont change current folder
parallel_batch_num = 2;
%10 is good for INTEL I-9 10850k(very hot), 30 is good for XEON w2295

for i = 1:1000
    folder_name = ['GA_fort' sprintf('%d',i)];
    if isfolder(folder_name)
       delete([pwd '\' folder_name '\*'])
       rmdir(folder_name)
       
    end
end
fname = 'start_run.bat';
fid = fopen(fname,'w');
for i = 1:parallel_batch_num
    folder_name = ['GA_fort' sprintf('%d',i)];
    mkdir(folder_name)
    ga_loc = [pwd '\GA_fort\run_ga.exe'];
    rpt_loc = [pwd '\input_data\*'];
    run_loc = [pwd '\' folder_name];
    copyfile(ga_loc,folder_name)
    copyfile(rpt_loc,folder_name)
    fprintf(fid,['start /d "',insertAfter(run_loc,'\','\') ,'" run_ga.exe\n']);
end
fclose('all')
system('start_run.bat')


%% 1.3 Shut down and collect data 

%you dont have to wait for independent population number to reach start_pt 
%you can just simply turn off the exe. and harvest
%To harvest: run this section 

%define your output name here, it will be a rec file (actually a txt but
%named .rec to differ from the other things), make sure you change it in
%each time of different run
%output name better record brief reason of run, although run info are
%written using fortran 
output_name = '79i_NA_3';
if isfolder('record')
else
    mkdir('record')
end
    

fname = 'collect_results.bat';
fid = fopen(fname,'w');
fprintf(fid,'@echo off \n');
str1 = ['set TxtPath=',insertAfter(pwd,'\','\'),'\\ \n'];
fprintf(fid, str1);
str2 = ['for /r %%TxtPath%% %%%%a in (*_out.txt) do more "%%%%a">> ' [insertAfter(pwd,'\','\'),'\\record\\'] output_name '.rec'];
fprintf(fid, str2);
fclose('all');    
system('collect_results.bat')


%% 1.3(alter) (run only if 1.3 failed) Shut down and collect data 

%you dont have to wait for independent population number to reach start_pt 
%you can just simply turn off the exe. and harvest
%To harvest: run this section 

%define your output name here, it will be a rec file (actually a txt but
%named .rec to differ from the other things), make sure you change it in
%each time of different run
%output name better record brief reason of run, although run info are
%written using fortran 
output_name = '79i_NA_6';
if isfolder('record')
else
    mkdir('record')
end
    

fname = 'collect_results2.bat';
fid = fopen(fname,'w');
fprintf(fid,'@echo off \n');
str1 = ['for /r "." %%%%a in (*_out.txt) do ( \n'];
fprintf(fid, str1);
str2 = ['type %%%%a >>'  [insertAfter(pwd,'\','\'),'\\record\\'] output_name '.rec \n'];
fprintf(fid, str2);
str3 = [') \n'];
fprintf(fid, str3);
fclose('all');    
system('collect_results2.bat')

%% 1.4 How to read the data 
%%Now you have the .rec file merging all the txt file stored in record
%folder, open it with Notepad, search in the note pad using the keyword
%[best obj(theta) popwise], try to find the least value, and the best material parameter set theta
%is written near the least value of popwise best obj(theta), EZ to do
%manually within 1 minute


%% 2-----------------start here if you choose to run it on cluster, linux based----------------------------------------------
%% 2.1 run this section to collect files needed to run on cluster 
if isfolder('cluster_run')
else
    mkdir('cluster_run')
end
delete([pwd '\cluster_run\*'])
rpt_loc = [pwd '\input_data\*'];
fortran_loc = [pwd '\GA_fort\*.f90'];
copyfile(rpt_loc,'cluster_run')
copyfile(fortran_loc,'cluster_run')
fclose('all')

%% 2.2
%1. select a base directory on cluster

%2. make a sub directory called GA_fort, copy things from cluster_run
%(local) directory onto GA_fort (cluster) directory

%3. copy bash files from cluster_bash (local) directory onto the base
%directory you chose on cluster 

%4. Open something like putty, run the bash file sequentially 
% from ga_0th.sh to ga_5th.sh

%5. When you want to terminate and harvest, either use scancel command,
%or use ga_6th.sh and modify the user name to your own 

%6. create an empty txt (with a space typed in) file named empty.txt in base folder

%7. run ga_7th.sh, then modify the yourchoiceofname.rec to the filename
%your want.rec and run ga_8th.sh  This gives you same things as running
%locally 




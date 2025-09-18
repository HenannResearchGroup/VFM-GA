% -------------VFM GA hyperelastic V1.0----------
% -------step1--1st func Data generation---------
% ------------------Zicheng Yan------------------
% -----Henann Research Group, Brown University---

%% First functionality, fit stress-strain and strain strain curve
%Given stress-strain and lateral-axial strain curve
%Generate input data for GA Hyper v1.0 fortran code package
%Sampling details based on description in GA hyper paper

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%look carefully into INPUT DATA section and Control Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
%% Input DATA
filename = 'sample_data.mat'; %change it to your filename
fakefile = 'poron_lo_quasi_data.mat'; %keep it

%%%Data structure requirement: compression and tension fully separated
%%%Six vectors in total(in double format), just the loading part
%%%1. compression_ax_stress (in Pa, engineering stress <=0 stress values)
%%%2. tension_ax_stress (in Pa, engineering stress >=0 stress values)
%%%3. compression_ax_strain (axial engineering strain)
%%%4. tension_ax_strain (axial engineering strain)
%%%5. compression_lat_strain(lateral engineering strain)
%%%6. tension_lat_strain (lateral engineering strain)
%%%1,3,5 should have same length, same goes for 2,4,6

input_data_mode = 2;
% =1 for complete data structure, (you have variable 1-6)
% =2 without tension, only compression ss and lateral curves,(you have variable 1,3,5)
% =3 only compression and tension ss curve, no lateral,(you have variable 1,2,3,4)
% =4 Only compression stress-strain curve, not even lateral, OOPS (you have variable 1,3)
% A.load your .mat file containing stress strain curve and lateral-axial strain
%curve

if input_data_mode == 1 
    %fantastic
    load(filename)
elseif input_data_mode == 2
    %due to the data structure, we add fake tension data from author's data
    %set, but we wont actually use tension in fortran GA hyper, Very LEGIT
    load(fakefile)
    clear compression_ax_stress
    clear compression_ax_strain
    clear compression_lat_strain
    load(filename)
elseif input_data_mode == 3
    %Fake up a lateral data, will address the fake part in fortran, not
    %that of a good choice, just go and get your DIC setup, kinda fine
    load(filename)
    clear compression_lat_strain
    clear tension_lat_strain
    max_lat_mag = 0.1; %we get a linear lat strain curve, change max mag to what you want(in compression)
    for i = 1:length(compression_ax_strain)
        compression_lat_strain(i) = linear_inter(0,0,compression_ax_strain(end),max_lat_mag,compression_ax_strain(i));
    end
    max_ten_lat_mag = max_lat_mag/abs(compression_ax_strain(end))*abs(tension_ax_strain(end));
    for i = 1:length(tension_ax_strain)
        tension_lat_strain(i) = -linear_inter(0,0,tension_ax_strain(end),max_ten_lat_mag,tension_ax_strain(i));
    end    
elseif input_data_mode == 4
    %Fake up nearly everything, youd better have better experimental data,
    %otherwise what you get out of package is not that physically
    %meaningful
    load(fakefile)
    clear compression_ax_stress
    clear compression_ax_strain
    clear compression_lat_strain
    load(filename)
    max_lat_mag = 0.1; %we get a linear lat strain curve, change max mag to what you want(in compression)
    for i = 1:length(compression_ax_strain)
        compression_lat_strain(i) = linear_inter(0,0,compression_ax_strain(end),max_lat_mag,compression_ax_strain(i));
    end
    max_ten_lat_mag = max_lat_mag/abs(compression_ax_strain(end))*abs(tension_ax_strain(end));
    for i = 1:length(tension_ax_strain)
        tension_lat_strain(i) = -linear_inter(0,0,tension_ax_strain(end),max_ten_lat_mag,tension_ax_strain(i));
    end    
end

%% CONTROL PARAMETERS
% B.number of samplings
%number of samplings in compression 
ncomp = 100; %(not recommended to change)
%number of samplings in tension
nten = 40; 
%(if first run use 35, for future run, if fitting in tension not desirable, increase the number, or vice versa)

%%%For fitting hyperelastic model:https://github.com/HenannResearchGroup/ElastomericFoam
%%%Because it is a model that focus on capturing the complex compression
%%%based material, when ncomp = 100, the choice of nten is usually
%%%between 25-45 recommended to be between 30-40 so fitting has a natural
%%%weight on compression. If you are trying to fit your own hyperelastic
%%%model, you can adjust the two numbers freely based on your own emphasis on
%%%the mechanical behavior in either tension/compression

% C.data format 
nel = 4; %(not recommended to change)
%according to empirical finding,9 nodes 4 elements work the best for first
%functionality of GA hyperelastic code package


% D.parameters to trim off tips of parameters
tstart = 4;
%number of data points to trim off from tension start(starting from origin)
cstart = 1;
%number of data points to trim off from compression start (starting from origin)

%%%You can trim off 3-5% eng strain in compression, bc G0 and K in 1st
%%%functionality need to be estimated first, trimming off doesn't affect the
%%%performance of GA hyper, sometimes experimental data at small strain using
%%%DIC could be a little messy

ttrim = 0;
%number of data points to trim off from tension end
ctrim = 0;
%number of data points to trim off from compression end

%%%Sometimes stress-strain curve could contain little bit of unloading and
%%%noise and the end tips at large deformation, ctrim and ttrim trims the
%%%curve giving you a nice thorough loading stress-strain and lat-axial
%%%strain curve
%% De Process

figure(2)
hold on 
plot(compression_ax_strain,compression_ax_stress,'-k','linewidth',2)
hold on 
plot(tension_ax_strain,tension_ax_stress,'-k','linewidth',2)
hold on

figure(1)
hold on
plot(tension_ax_strain,tension_lat_strain,'-k','linewidth',2)
hold on 
plot(compression_ax_strain,compression_lat_strain,'-k','linewidth',2)
hold on

    
[mins,mindex] = max(abs(compression_ax_strain));
[maxs,maxdex] = max(tension_ax_stress);



mindex = mindex-ctrim;
maxdex = maxdex-ttrim;



ncomp_stress = compression_ax_stress./compression_ax_stress(mindex);
ncomp_strain = compression_ax_strain./compression_ax_strain(mindex);

nten_stress = tension_ax_stress./tension_ax_stress(maxdex);
nten_strain = tension_ax_strain./tension_ax_strain(maxdex);
samples = ncomp+nten;
comp_dist = zeros(1,mindex+1-cstart);
acc_comp_dist = zeros(1,mindex+1-cstart);
for i = cstart+1:mindex
    index = i-cstart+1;
    comp_dist(index) = sqrt((ncomp_stress(i)-ncomp_stress(i-1))^2+(ncomp_strain(i)-ncomp_strain(i-1))^2);
end

for i = 1:length(comp_dist)
    start = 1;
    while start<=i
        acc_comp_dist(i) = acc_comp_dist(i)+comp_dist(start);
        start = start+1;
    end
end

dist_inc = acc_comp_dist(end)/(ncomp);
for i = 1:ncomp
    dist = dist_inc*(i-1);
    start = 1;
    while dist>acc_comp_dist(start)
        start = start+1;
    end
    if start-1 <= 0 && dist == 0
        index_out(i) = start;
    else
        index_out(i) = linear_inter(acc_comp_dist(start-1),start-1,acc_comp_dist(start),start,dist);
    end
end

index_out = index_out+cstart-1;




ten_dist = zeros(1,maxdex+1-tstart);
acc_ten_dist = zeros(1,maxdex+1-tstart);
for i = tstart+1:maxdex
    index = i-tstart+1;
    ten_dist(index) = sqrt((nten_stress(i)-nten_stress(i-1))^2+(nten_strain(i)-nten_strain(i-1))^2);
end

for i = 1:length(ten_dist)
    start = 1;
    while start<=i
        acc_ten_dist(i) = acc_ten_dist(i)+ten_dist(start);
        start = start+1;
    end
end

ten_dist_inc = acc_ten_dist(end)/(nten);
for i = 1:nten
    dist = ten_dist_inc*(i-1);
    start = 1;
    while dist>acc_ten_dist(start)
        start = start+1;
    end
    if start-1 <= 0 && dist == 0
        tindex_out(i) = start;
    else
        tindex_out(i) = linear_inter(acc_ten_dist(start-1),start-1,acc_ten_dist(start),start,dist);
    end
end

tindex_out = tindex_out+tstart-1;



    
    
if nel == 4
    coord(:,1) = [1,2,3,4,5,6,7,8,9]';
    coord(:,2) = [0,0.5,1,0,0.5,1,0,0.5,1]';
    coord(:,3) = [0,0,0,0.5,0.5,0.5,1,1,1]';
    element = [1,1,2,5,4;2,2,3,6,5;3,4,5,8,7;4,5,6,9,8];
   fileID1 = strcat('fit9_node.inp');
    fileID2 = strcat('fit9_element.inp');
    format1 = '%d , %8.6f , %8.6f\n';
    format2 = '%d, %d, %d, %d, %d\n';    
    fileprint1 = fopen(fileID1,'w');
    fileprint2 = fopen(fileID2,'w');
    for i = 1:9
        vec = [coord(i,1),coord(i,2),coord(i,3)];
        fprintf(fileprint1,format1,vec);
    end
    for i = 1:4
        fprintf(fileprint2,format2,element(i,:));
    end
    
    for i = 1:ncomp
        index = index_out(i);
        if index == round(index)
            as = compression_ax_strain(index);
            ls = compression_lat_strain(index);
            astress = compression_ax_stress(index);

        else
            fi = floor(index);
            ci = ceil(index);
            as1 = compression_ax_strain(fi);
            as2 = compression_ax_strain(ci);
            as = linear_inter(fi,as1,ci,as2,index);
            ls1 = compression_lat_strain(fi);
            ls2 = compression_lat_strain(ci);
            ls = linear_inter(fi,ls1,ci,ls2,index);
            astress1 = compression_ax_stress(fi);
            astress2 = compression_ax_stress(ci);
            astress = linear_inter(fi,astress1,ci,astress2,index);

        end
        ydisp = as;
        xdisp = ls;
        force = astress;
        rpt(1,:) = [1,0,0,-astress/3,0,0,0];
        rpt(2,:) = [2,0,0,-astress/3,0,xdisp/2,0];
        rpt(3,:) = [3,0,0,-astress/3,0,xdisp,0];
        rpt(4,:) = [4,0,0,0,0,0,ydisp/2];
        rpt(5,:) = [5,0,0,0,0,xdisp/2,ydisp/2];
        rpt(6,:) = [6,0,0,0,0,xdisp,ydisp/2];
        rpt(7,:) = [7,0,0,astress/3,0,0,ydisp];
        rpt(8,:) = [8,0,0,astress/3,0,xdisp/2,ydisp];
        rpt(9,:) = [9,0,0,astress/3,0,xdisp,ydisp];
        format = '%d     %8.8f     %8.8f     %8.8f     %8.8f     %8.8f     %8.8f\n';
        if samples<=100
            if i<11
                fileID = strcat('w','0',int2str(i-1),'.rpt');
            else
                fileID = strcat('w',int2str(i-1),'.rpt');
            end
        else
            if i<11
                fileID = strcat('w','00',int2str(i-1),'.rpt');
            elseif i<101
                fileID = strcat('w0',int2str(i-1),'.rpt');
            else
                fileID = strcat('w',int2str(i-1),'.rpt');
            end
        end
        fileprint = fopen(fileID,'w');
        for j = 1:9
            fprintf(fileprint,format,rpt(j,:));
        end
        figure(2)
        hold on
        plot(as,astress,'r*','LineWidth',3,'markersize',8)
        hold on
        
        figure(1)
        hold on
        plot(as,ls,'r*','LineWidth',3,'markersize',8)
        hold on
    end
    start = ncomp;
    for i = 1:nten
        index = tindex_out(i);
        if index == round(index)
            as = tension_ax_strain(index);
            ls = tension_lat_strain(index);
            astress = tension_ax_stress(index);
        else
            fi = floor(index);
            ci = ceil(index);
            as1 = tension_ax_strain(fi);
            as2 = tension_ax_strain(ci);
            as = linear_inter(fi,as1,ci,as2,index);
            ls1 = tension_lat_strain(fi);
            ls2 = tension_lat_strain(ci);
            ls = linear_inter(fi,ls1,ci,ls2,index);
            astress1 = tension_ax_stress(fi);
            astress2 = tension_ax_stress(ci);
            astress = linear_inter(fi,astress1,ci,astress2,index);
        end
        ydisp = as;
        xdisp = ls;
        force = astress;
        rpt(1,:) = [1,0,0,-astress/3,0,0,0];
        rpt(2,:) = [2,0,0,-astress/3,0,xdisp/2,0];
        rpt(3,:) = [3,0,0,-astress/3,0,xdisp,0];
        rpt(4,:) = [4,0,0,0,0,0,ydisp/2];
        rpt(5,:) = [5,0,0,0,0,xdisp/2,ydisp/2];
        rpt(6,:) = [6,0,0,0,0,xdisp,ydisp/2];
        rpt(7,:) = [7,0,0,astress/3,0,0,ydisp];
        rpt(8,:) = [8,0,0,astress/3,0,xdisp/2,ydisp];
        rpt(9,:) = [9,0,0,astress/3,0,xdisp,ydisp];
        format = '%d     %8.8f     %8.8f     %8.8f     %8.8f     %8.8f     %8.8f\n';
        if samples<=100
            if i+start<11
                fileID = strcat('w','0',int2str(i+start-1),'.rpt'); 
            else
                fileID = strcat('w',int2str(i+start-1),'.rpt');
            end
        else
            if i+start<11
                fileID = strcat('w','00',int2str(i+start-1),'.rpt');
            elseif i+start<101
                fileID = strcat('w0',int2str(i+start-1),'.rpt');
            else
                fileID = strcat('w',int2str(i+start-1),'.rpt');
            end
        end
        fileprint = fopen(fileID,'w');
        for j = 1:9
            fprintf(fileprint,format,rpt(j,:));
        end
        figure(2)
        hold on 
        plot(as,astress,'r*','LineWidth',3,'markersize',8)
        hold on
        
        figure(1)
        hold on 
        plot(as,ls,'r*','LineWidth',3,'markersize',8)
        hold on
    end
end

fileID3 = strcat('record_input_data_setup.txt');
format3 = '%d   %d   %d\n';
fileprint3 = fopen(fileID3,'w');
fprintf(fileprint3,format3,[input_data_mode ncomp nten]);

fileID4 = strcat('record_expdata_name.txt');
fileprint4 = fopen(fileID4,'w');
fprintf(fileprint4,filename);


delete([pwd '\input_data\*'])
copyfile('*.rpt', [pwd '\input_data'])
fclose('all');
delete([pwd '\*.rpt'])


    
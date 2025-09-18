% ---------------VFM GA hyperelastic V1.0---------------
% -------2nd functionality input data structure---------
% -------------------Zicheng Yan------------------------
% -----Henann Research Group, Brown University----------



%% Second functionality, fit processed field data from DIC experiment
%Given an undeformed configuration (mesh) and deformed configurations 
%The 2nd functionality fits the DIC field data for constitutive model


%%%%%%%%%%%%%%%%%%%%%%%%
%%%This is not a code, but an instruction for data structure of DIC
%%%experiment, prepare the required files and put into /input_data folder
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%In the input_data folder, an example data for the 2nd funtionality first
%%%example in the paper is stored and you can play around with it.
%%%In S_1 matlab file, the input parameters are setup for the example data.
%%%You can use the example to understand the basics of the code.

%% 
%1. the fortran code allow utilizing results from 2 DIC experiment together

% For the first configuration (experiment)
% You need comp_node.inp, comp_element.inp as the node and element file 
% To define the initial mesh by DIC

% 1.1 Node file data structure (same as abaqus node part in inp file)
% node_number, position_in_x, position_in_y
% ........
% ........


% 1.2 Element file data structure (same as abaqus element part in inp file)
% element_number, node_number1, node_number2, node_number3, node_number4
% ........
% ........
%(4 noded square element strongly recommended) 


% you need to define deformed config using file named after w recording the
% displacement measured by DIC and reaction forces measured by load cell,
% starting with w00.rpt with numbering from small to large absolute nominal strain. 


% 1.3 rpt file data structure (same as the field output file of abaqus when
% asking for reaction force and displacement) 
% node_number,0,RF in x, RF in y, 0, disp in x, disp in y
% .............
% .............
%(Note RF refers to reaction force, disp refers to displacement) 
%(Since reaction force is measured as a scalar in normal Instron
%compression/tension experiment, to convert measured force scalar to put
%into the rpt files, simply divide the value by the number of nodes on the
%top or bottom and assign that value to top/bottom nodes only, see the
%example w00.rpt file as an example)


%2. for the second configuration, you need ten_node.inp and ten_element.inp
%to define the undeformed configuration, and t00.rpt, t01.rpt, etc to
%define the deformed configurations. Data structure stays the same. 

%3. If you are fitting using near simple compression-tension field data with 2 configurations
%Fill comp_node.inp, comp_element.inp, w00.rpt - w??.rpt using compression
%data, and ten_node.inp, ten_element.inp, t00.rpt - t??.rpt using tension
%data. 

% If you are fitting only using 1 configuration (fit to inhomogeneous field data or half response only)
% You can just have comp_node.inp comp_element.inp and the w00.rpt -
% w??.rpt

%4. To fit to more than 3 experiments with field output together, you need
%to understand and modify the code, stack allocation strongly recommended,
%dont try to convert input preprocessor function to heap, you would regret
%it due to the loss in speed

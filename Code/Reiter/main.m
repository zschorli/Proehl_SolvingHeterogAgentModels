for case_nr=1:4
global Params;
Params=[];
Params.case_nr = case_nr;    
delete('res/*')

% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% main program:
jedcpath = pwd;

% ONLY NECESSARY ONCE AFTER INSTALLATION:
% mkdir res;
makeshocks;
% GENERATE DLL-FILES (ONLY NECESSARY IF INCLUDED DLL-s DON'T WORK)
% cd libm; allmex; cd(jedcpath);

% SET PATHS:
addpath([jedcpath '/libm/']);
addpath([jedcpath '/lib_mirandafackler/']);
addpath([jedcpath '/data/']);

comp_time = tic;
% SOLVE FOR ALL PARAMETER COMBINATIONS:
allsolve;
comp_time = toc(comp_time);

fn = resfilename('V',maketag(Params),1);
load(fn,'VAac','xEqu','Params2');
save(strcat('R_Sol',num2str(case_nr),'.mat'),'VAac','xEqu','Params2','comp_time')

% % SIMULATE FOR ALL PARAMETER COMBINATIONS (GENERATES OUTPUT NECESSARY FOR COMPARISON):
% allsimul;
% 
% % GENERATE TABLE 1 OF PAPER:
% tablesimul;


end
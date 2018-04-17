function ea_loaddti_exploredti(varargin)
% ea_loaddti_exploredti(inputdir), ea_loaddti_exploredti()
% 
% user is asked to save the *.mat TRACT file and the *.nii T1 file from the
% processed ExploreDTI files. The selected filenames are then saved to the
% ea_exploredti_file.mat file.

% Function created by Wouter Potters, Academic Medical Center, Amsterdam

if nargin > 0 
    inputdir = varargin{1};
else
    inputdir = cd;
end

[file,directory] = uigetfile(fullfile(inputdir,'*.mat'),'Select *.mat TRACT file');
dti_file = fullfile(directory,file); 

[file,directory] = uigetfile(fullfile(directory,'*.nii'),'Select *.nii T1 file belonging to the DTI pipeline.');
t1_file = fullfile(directory,file); 

save(fullfile(inputdir,'ea_exploredti_file.mat'),'dti_file','t1_file'); %#ok<*NASGU>

end
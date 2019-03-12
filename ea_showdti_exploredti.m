function ea_showdti_exploredti(resultfig,options)
% - loads the ea_exploredti_file.mat file
% - performs a registration between t1 from exploredti and t1 in leaddbs
% - transforms the tracts calculated in exploredti to the leaddbs coordinate 
%   space.
% - save these distances to *_distances.txt
% - plot histogram for each distance distribution


% Function created by Wouter Potters, Academic Medical Center, Amsterdam


options = getappdata(resultfig,'options');
patdir = options.uipatdirs;

exploredti_settings_file = fullfile(patdir{1},'ea_exploredti_file.mat');
exploredti_settings = load(exploredti_settings_file);

% hold on;
tract = load(exploredti_settings.dti_file);
% nii_t1_exdti = nifti(exploredti_settings.t1_file);
% nii_t1_leaddbs = nifti(fullfile(options.uipatdirs{1}, options.prefs.prenii_unnormalized_t1));
tract.Tracts_lead = cell(1,length(tract.Tracts));

ref = fullfile(options.uipatdirs{1}, options.prefs.prenii_unnormalized_t1);
moving = [exploredti_settings.t1_file(1:end-4) '_leaddbsREG.nii'];
if exist([exploredti_settings.t1_file(1:end-4) '_leaddbsREG.nii'],'file') ~= 2
    try
        copyfile(exploredti_settings.t1_file,moving);
    catch err
        error('Copying file %s failed',[exploredti_settings.t1_file(1:end-4) '_leaddbsREG.nii'])
    end
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {ref};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {moving};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'ncc';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);
end
% orig_mat = nii_t1_exdti.mat; % original T1 nii transformation matrix
moved_mat = nifti(moving); % explore dti transformation matrix
moved_mat = moved_mat.mat; % explore dti transformation matrix registered to LEAD DBS

for ind = 1:length(tract.Tracts)
    % step 1: convert exploredti xyz to ijk
    tmat = tract.VDims;
    tijk = tract.Tracts{ind}./tmat;
    tijk(:,1) = size(tract.TractMask,1) - tijk(:,1);
    tijk(:,2) = size(tract.TractMask,2) - tijk(:,2);
    tijk1 = [tijk ones(size(tijk,1),1)];
    txyz1 = ( ([0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1]) * tijk1')'; % scale to include plot scale of 2x
    txyz1 = ((moved_mat) * txyz1')';
    tract.Tracts_lead{ind} = txyz1(:,[1 2 3]);
end

% EXPLOREDTI:
% mm coordinates for tracts
% flip lr J! (radiological convention)
% jj = max(j) - j + 1
% xyz = TM * i(jj)k'

% LEAD_DBS:
% ijk = 0.5:size+0.5
% xyz' = TM * ijk'


%%% PLOT DATA
hold on;
all_tracts = cellfun(@(x) patch('XData',[x(:,1); nan],...
    'YData',[x(:,2); nan],...
    'ZData',[x(:,3); nan] ,...
    'FaceVertexCdata',x(:,3),...
    'FaceColor','none','edgecolor','flat'),...
    tract.Tracts_lead,'uniformoutput',0);
colormap(jet(100)); c = colorbar; set(c,'color','w'); t= title(c,'Distance to contactpoint (mm)'); set(t,'color','w')
caxis([0 100]) % in mm

cf = get(gcf,'Number');
try % close all windows that are not the Electrode scene.
    all_fig = cell2mat(get(findobj('Name','DICOM: Electrode-Scene...building...'),'Number')); % find electrode scene
    close(all_fig(~ismember(all_fig,cf))) % close all windows that are not the Electrode scene.
end

active_window = findobj('Name','DICOM: Electrode-Scene...building...');
opts = get(active_window,'Children');
active_axis = opts(arrayfun(@(x) isa(x,'matlab.graphics.axis.Axes'),opts));
allplots = get(active_axis,'Children');
patches = allplots(arrayfun(@(x) isa(x,'matlab.graphics.primitive.Patch'),allplots));
% arrayfun(@(x) get(x),patches) % not used


%%% CALCULATE DISTANCE TO ELECTRODES
load(fullfile(patdir{1},'ea_reconstruction.mat'))
fid = fopen(fullfile([exploredti_settings.dti_file(1:end-4) '_distances.txt']),'w');
for selected_electrode = 1:2
    fprintf(fid,'\n\n====================================');
    fprintf(fid,'\nelectrode %i',selected_electrode);
    for point = 1:4
        fprintf(fid,'\npoint %i',point);
        fprintf(fid,'\n%f %f %f',reco.native.coords_mm{selected_electrode}(point,:));
        for ind = 1:length(all_tracts)
            dst = sqrt(sum((reco.native.coords_mm{selected_electrode}(point,:) - all_tracts{ind}.Vertices).^2,2));
            all_tracts{ind}.FaceVertexCData = dst;
        end
        dst_all = cellfun(@(x) x.FaceVertexCData,all_tracts,'uniformoutput',false);
        
        % only electrodes that are close (< 1 cm)
        dst_all_selection = cellfun(@(x) any(x.FaceVertexCData < 10),all_tracts,'uniformoutput',true);
        dst_shortest = cellfun(@(x) min(x), dst_all(dst_all_selection));
        fprintf(fid,'\nnumber of fibers: %i',numel(dst_shortest));
        fprintf(fid,'\nminimal and maximal (10mm max) distance to nearby (<10 mm) fibers: %f to %f',min(dst_shortest),max(dst_shortest));
        fprintf(fid,'\nmedian +/- SD distance to nearby (<10 mm) fibers: %f +/- %f',median(dst_shortest),std(dst_shortest));
        
        %%% PLOT DISTANCE TO ELECTRODES IN HISTOGRAM
        fg=gcf;
        figure('Name',sprintf('electrode %i, point %i',selected_electrode,point)); hist(dst_shortest);
        title(sprintf('Histogram of distance to nearby fibers (< 10 mm), #%i',numel(dst_shortest)));
        figure(fg); % go back to leaddbs figure
        
        fprintf('\n')
        fprintf('\n---')
    end
end
fclose(fid);
function ea_dti_show_exploredti(resultfig)
% - loads the ea_exploredti_file.mat file
% - performs a registration between t1 from exploredti and t1 in leaddbs
% - transforms the tracts calculated in exploredti to the leaddbs coordinate 
%   space.
% - save these distances to *_distances.txt
% - plot histogram for each distance distribution


% Function created by Wouter Potters, Academic Medical Center, Amsterdam

MAXIMAL_DISTANCE_TO_INCLUDE_IN_HISTOGRAMS = 15; % mm

options = getappdata(resultfig,'options');
patdir = options.uipatdirs;

if options.atl.ptnative ~= 1
    delete(resultfig)
    error('DTI visualisation only works in native patient space. Select Native Space in item 8, first dropdown.')
end

exploredti_settings_file = fullfile(patdir{1},'ea_exploredti_file.mat');
exploredti_settings = load(exploredti_settings_file);

% hold on;
tract = load(exploredti_settings.dti_file);
% nii_t1_exdti = nifti(exploredti_settings.t1_file);
% nii_t1_leaddbs = nifti(fullfile(options.uipatdirs{1}, options.prefs.prenii_unnormalized_t1));
tract.Tracts_lead = cell(1,length(tract.Tracts));

dti = fullfile([exploredti_settings.t1_file(1:end-10) '.nii']);
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
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);
end

% orig_mat = nii_t1_exdti.mat; % original T1 nii transformation matrix, not used
moved_mat = nifti(moving); % explore dti transformation matrix
moved_mat = moved_mat.mat; % explore dti transformation matrix registered to LEAD DBS
for ind = 1:length(tract.Tracts)
    % step 1: convert exploredti xyz to ijk
    tmat = tract.VDims;
    tijk = tract.Tracts{ind}./tmat;
    tijk(:,1) = size(tract.TractMask,1) - tijk(:,1);
    tijk1 = [tijk ones(size(tijk,1),1)];
    txyz1 = ( ([0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1]) * tijk1')';
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

%%% CALCULATE DISTANCE TO ELECTRODES
load(fullfile(patdir{1},'ea_reconstruction.mat'))
fid = fopen(fullfile([exploredti_settings.dti_file(1:end-4) '_distances.txt']),'w');
for selected_electrode = 1:2
    fprintf(fid,'\n\n====================================');
    fprintf(fid,'\nelectrode %i',selected_electrode);
    for point = 1:4
        fprintf(fid,'\npoint %i',point);
        fprintf(fid,'\n%f %f %f',reco.native.coords_mm{selected_electrode}(point,:));
        for ind = 1:length(tract.Tracts_lead)
            dst = sqrt(sum((reco.native.coords_mm{selected_electrode}(point,:) - tract.Tracts_lead{ind}).^2,2));
            tract.Tracts_lead_dst{ind} = dst;
        end
        dst_all = cellfun(@(x) x,tract.Tracts_lead_dst,'uniformoutput',false);
        
        % only electrodes that are close (< 1 cm)
        dst_all_selection = cellfun(@(x) any(x < MAXIMAL_DISTANCE_TO_INCLUDE_IN_HISTOGRAMS),tract.Tracts_lead_dst,'uniformoutput',true);
        dst_shortest = cellfun(@(x) min(x), dst_all(dst_all_selection));
        fprintf(fid,'\nnumber of fibers: %i',numel(dst_shortest));
        fprintf(fid,'\nminimal and maximal (10mm max) distance to nearby (<10 mm) fibers: %f to %f',min(dst_shortest),max(dst_shortest));
        fprintf(fid,'\nmedian +/- SD distance to nearby (<10 mm) fibers: %f +/- %f',median(dst_shortest),std(dst_shortest));
        
        %%% PLOT DISTANCE TO ELECTRODES IN HISTOGRAM
        fg=gcf;
        f1 = figure('Name',sprintf('electrode %i, point %i',selected_electrode,point)); hist(dst_shortest);
        set(f1,'Tag','electrode_histogram');
        title(sprintf('Histogram of distance to nearby fibers (< %1.1f mm), #%i',MAXIMAL_DISTANCE_TO_INCLUDE_IN_HISTOGRAMS,numel(dst_shortest)));
        axxxx = axis; axis(gca, [0 MAXIMAL_DISTANCE_TO_INCLUDE_IN_HISTOGRAMS 0 axxxx(4)]);
        figure(fg); % go back to leaddbs figure
        
        fprintf('\n')
        fprintf('\n---')
        if point == 1
            dst_all_pts = dst_all;
        else
            dst_all_pts = cellfun(@(x,y) min([x,y],[],2),dst_all,dst_all_pts,'uniformoutput',false);
        end
    end
    if selected_electrode == 1
        dst_all_pts_el = dst_all_pts;
    else
        dst_all_pts_el = cellfun(@(x,y) min([x,y],[],2),dst_all_pts,dst_all_pts_el,'uniformoutput',false);
    end
end
tract.Tract_lead_dst = dst_all_pts_el;
fclose(fid);

%%% PLOT DATA
hold(gca(resultfig),'on');
all_tracts = cellfun(@(x,y) patch(gca(resultfig),'XData',[x(:,1); nan],...
    'YData',[x(:,2); nan],...
    'ZData',[x(:,3); nan] ,...
    'FaceColor','flat','EdgeColor','interp','facevertexcdata',[y; nan]),...
    tract.Tracts_lead,tract.Tract_lead_dst,'uniformoutput',0);
colormap(gca(resultfig),[jet(100); 1 1 1]); c = colorbar; set(c,'color','w'); t= title(c,{'Minimal distance to';'contactpoint (mm)'; sprintf('(max. %1.0f mm)',MAXIMAL_DISTANCE_TO_INCLUDE_IN_HISTOGRAMS)}); set(t,'color','w')
caxis(gca(resultfig),[0 MAXIMAL_DISTANCE_TO_INCLUDE_IN_HISTOGRAMS]) % in mm

cf = get(gcf,'Number');
try % close all windows that are not the Electrode scene.
    all_fig = cell2mat(get(findobj('Name','DICOM: Electrode-Scene...building...'),'Number')); % find electrode scene
    close(all_fig(~ismember(all_fig,cf))) % close all windows that are not the Electrode scene.
end

% active_window = findobj('Name','DICOM: Electrode-Scene...building...');
% opts = get(active_window,'Children');
% active_axis = opts(arrayfun(@(x) isa(x,'matlab.graphics.axis.Axes'),opts));
% allplots = get(active_axis,'Children');
% patches = allplots(arrayfun(@(x) isa(x,'matlab.graphics.primitive.Patch'),allplots));
% arrayfun(@(x) get(x),patches) % not used
% BATCH REGISTER ALL DATA
check = true;

ldbs_dir = 'L:\basic\divd\knf\Onderzoek_studenten\LauraKoster\5. electrode reconstruction pacer default (9, 21, 25, 29 zonder bsc)';
expdti_dir = 'L:\basic\divd\knf\Onderzoek_studenten\LauraKoster\nieuwe DTI data';
lead_dbs_dir = dir(ldbs_dir);
lead_dbs_dir = lead_dbs_dir(3:end);
t1names_leaddbs = cellfun(@(x) fullfile(ldbs_dir,x,'anat_t1.nii'),{lead_dbs_dir([lead_dbs_dir.isdir]).name},'uniformoutput',false);
t1folders = cellfun(@(x) fullfile(expdti_dir,regexp(x,'PD_DBS[0-9]+','match','once'),[regexp(x,'PD_DBS[0-9]+','match','once') '_DTI'],'*T1.nii'),{lead_dbs_dir([lead_dbs_dir.isdir]).name},'uniformoutput',false);

for iPat = 1:length(t1folders)
    ref = t1names_leaddbs{iPat};
    dti = strrep(t1folders{iPat},'*T1.nii',ls(t1folders{iPat}));
    fprintf('Patient %i:\nexploredti: %s\nleaddbs: %s\n\n',iPat,dti,ref)
    moving = [exploredti_settings.t1_file(1:end-4) '_leaddbsREG.nii'];
    %try delete(moving); end % TRY TO DELETE ALL DATA FIRST
    if exist([exploredti_settings.t1_file(1:end-4) '_leaddbsREG.nii'],'file') ~= 2
        try
            copyfile(exploredti_settings.t1_file,moving);
        catch err
            error('Copying file %s failed',[exploredti_settings.t1_file(1:end-4) '_leaddbsREG.nii'])
        end
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {ref};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {moving};
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'mi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
    end
    
    if check
        spm_check_registration(ref,moving);
    end

end
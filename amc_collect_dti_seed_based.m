main_dir = 'L:\basic\divd\knf\Onderzoek_projecten\MT-DBS\Data\Lead DBS analyse met t1 en t2\6 VAT';
patients = dir(main_dir);
patients = patients(cellfun(@(p) ~isempty(regexp(p,'PD_DBS*','once')),{patients.name},'uniformoutput',true));

nifti_atlas = 'L:\basic\divd\knf\Onderzoek_projecten\MT-DBS\leadDBS_v2.2\templates\space\MNI_ICBM_2009b_NLIN_ASYM\labeling\Automated Anatomical Labeling (Tzourio-Mazoyer 2002).nii';
txt_atlas = [nifti_atlas(1:end-3) 'txt'];
[~,atlasname,~] = fileparts(nifti_atlas);
nifti_atlas = nifti(nifti_atlas);
atlas_labels = textscan(fopen(txt_atlas,'r'),'%f\t%s');
if ~all(diff(atlas_labels{1})==1), error('assumption that all labels are +1 increasing is invalid'); end
atlas_labels = atlas_labels{2}; % assuming that labels increase by 1 every time.
data_atlas = nifti_atlas.dat(:,:,:);

excel_export_file = ['results_VATDTI_atlas_' atlasname];
number2letter = @(n)char(n-1+'a');

% loop over all patient folders
for p = 1:length(patients)
    patient = patients(p).name;
    lead_dbs_stim_folder = fullfile(main_dir,patient,'stimulations');
    if ~isfolder(lead_dbs_stim_folder)
        warning('patient %s has no stimulations',patient)
        continue
    end
    
    % get all stimulation folders for current patient
    stimulations = dir(lead_dbs_stim_folder);
    stimulations = stimulations(cellfun(@(s) ~isempty(regexp(s,'^[LR]H_K.*mA$','once')),{stimulations.name}));
    for s = 1:length(stimulations)
        stimulation = stimulations(s).name;
        fprintf('processing patient %s, stimulation %s\n',patient, stimulation)
        nifti_results_tracts = dir(fullfile(lead_dbs_stim_folder,stimulation,'PPMI_90 (Ewert 2017)','*.nii'));
        if length(nifti_results_tracts) ~= 1
            error('size nifti_results_tracts: %i',length(nifti_results_tracts))
        end
        vat_left_struc_seed.nii
        nifti_tracts = nifti(fullfile(lead_dbs_stim_folder,stimulation,'PPMI_90 (Ewert 2017)',nifti_results_tracts.name));
        
        % calculate overlap for current stimulation in current patient
        sums = cell(length(atlas_labels),1);
        percentages = cell(length(atlas_labels),1);
        binaries = cell(length(atlas_labels),1);
        data_image = nifti_tracts.dat(:,:,:);
        for label = 1:length(atlas_labels)
            mask = (data_atlas == label);
            masked = (data_image .* mask);
            sums{label} = sum(masked(:));
            if isinf(sum(masked(:)~=0)/sum(mask(:))*100)
                percentages{label} = 0;
            else
                percentages{label} = sum(masked(:)~=0)/sum(mask(:))*100;
            end
            binaries{label} = any(masked(:) > 0);
        end
        
        % save results
        sheetname_excel = sprintf('%s', stimulation);
        warning off;
        
        % always write atlas labels
        xlswrite(fullfile(main_dir,[excel_export_file '_percent.xlsx']),[{''}; atlas_labels],sheetname_excel,'A1');
        xlswrite(fullfile(main_dir,[excel_export_file '_sums.xlsx']),[{''}; atlas_labels],sheetname_excel,'A1');
        xlswrite(fullfile(main_dir,[excel_export_file '_binaries.xlsx']),[{''}; atlas_labels],sheetname_excel,'A1');
        
        % always write patient numbers
        xlswrite(fullfile(main_dir,[excel_export_file '_percent.xlsx']),[{''} {patients.name}],sheetname_excel,'A1');
        xlswrite(fullfile(main_dir,[excel_export_file '_sums.xlsx']),[{''} {patients.name}],sheetname_excel,'A1');
        xlswrite(fullfile(main_dir,[excel_export_file '_binaries.xlsx']),[{''} {patients.name}],sheetname_excel,'A1');
        
        if p+1 > 26
            if p + 1 - 26 > 26
                error('fix me')
            end
            pn = ['a' number2letter(rem(p+1,26))];
        else
            pn = number2letter(p + 1);
        end
        xlswrite(fullfile(main_dir,[excel_export_file '_percent.xlsx']),[{patient}; percentages],sheetname_excel,sprintf('%s1',(pn)));
        xlswrite(fullfile(main_dir,[excel_export_file '_sums.xlsx']),[{patient}; sums],sheetname_excel,sprintf('%s1',(pn)));
        xlswrite(fullfile(main_dir,[excel_export_file '_binaries.xlsx']),[{patient}; binaries],sheetname_excel,sprintf('%s1',(pn)));
        warning on;
        fprintf('\b done\n')
    end
end
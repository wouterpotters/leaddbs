function ea_init_coregmrpopup(handles,refine)
if ~exist('refine','var')
    refine=0;
end

    cmethods={'SPM',...
        'FSL FLIRT',...
        'FSL BBR',...
        'ANTs',...
        'BRAINSFIT',...
        'Hybrid SPM & ANTs',...
        'Hybrid SPM & FSL',...
        'Hybrid SPM & BRAINSFIT'};
set(handles.coregmrpopup,'String',cmethods)
set(handles.coregmrpopup,'Value',1); % default SPM

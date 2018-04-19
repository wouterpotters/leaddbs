function h = ea_dti_addinterface(h)
% h = ea_dti_addinterface(h) where h is a handles structure.
%
% This function programmatically adds the DTI interface to Lead-DBS 
% main GUI defined by lead_dbs.fig. It should be called once (using opening
% function in lead_dbs.m)

% This function only exists to make GIT integration easy without changing 
% the existing *.fig file in the repository.
%
% author: Wouter Potters, Academic Medical Center, Amsterdam
% date: 17-4-2018

objects_have_been_moved = h.vizpanel.UserData; % only move the gui objects once.
narginchk(1,1); % requires only 1 input argument, 1 output argument.
nargoutchk(1,1); % requires only 1 input argument, no output arguments.

%% move objects to make space for DTI stuff.
objects_to_move = [h.run_button ...
                   h.exportcode ...
                   h.viewmanual ...
                   h.openpatientdir ...
                   h.vizpanel ...
                   h.overwriteapproved ...
                   h.statusone ...
                   h.statustwo ];

distance_move = 55;
               
if isempty(objects_have_been_moved) || ~objects_have_been_moved
    movedown(objects_to_move, -distance_move); % move by 20 pixels
    h.vizpanel.UserData = true;
end

%% rename visualization box from 7 to 8
h.vizpanel.Title(1) = '8';

%% add new panel for DTI
pos_psa = h.psapanel.Position;
pos_viz = h.vizpanel.Position;
pos_dti = [pos_psa(1) pos_viz(2)+pos_viz(4)+1 pos_psa(3) pos_psa(2)-(pos_viz(2)+pos_viz(4)+1)];
h.dtipanel = uipanel(h.leadfigure,...
    'Units','pixels',...
    'Position',pos_dti,...
    'Title','7. Load and add DTI fibers');
cellfun(@(c) set(h.dtipanel,c,get(h.vizpanel,c)),{'ForegroundColor','FontName','FontSize','FontUnits','BackgroundColor','ShadowColor'})

%% add DTI buttons in DTI panel with callback functions
h.chk_ExploreDTI_data = uicontrol(h.dtipanel,'Style','checkbox',...
    'String','ExploreDTI data',...
    'ButtonDownFcn',str2func('@(hObject,eventdata)lead_dbs(''include_exploredti_ButtonDownFcn'',hObject,eventdata,guidata(hObject))'),...
    'Callback',str2func('@(hObject,eventdata)lead_dbs(''include_exploredti_Callback'',hObject,eventdata,guidata(hObject))'));

h.btn_ExploreDTI_settings = uicontrol(h.dtipanel,'Style','pushbutton',...
    'String','Settings',...
    'Callback',str2func('@(hObject,eventdata)lead_dbs(''openexploredtidata_Callback'',hObject,eventdata,guidata(hObject))'));

cellfun(@(c) set(h.chk_ExploreDTI_data,c,get(h.normcheck,c)),{'ForegroundColor','FontName','FontSize','FontWeight','FontUnits','BackgroundColor','Position'})
cellfun(@(c) set(h.btn_ExploreDTI_settings,c,get(h.openleadconnectome,c)),{'ForegroundColor','FontName','FontWeight','FontSize','FontUnits','BackgroundColor','Position'})

function movedown(objects, pixel)
% this helper function moves objects down by 'pixel' pixels
narginchk(2,2);
for o = objects
    try % updating position for each object
        position = get(o,'Position');
        set(o,'Position',position + [0 pixel 0 0]);
    catch err
        warning('could not find object: %s (error: %s)', o.String, err.message);
    end
end
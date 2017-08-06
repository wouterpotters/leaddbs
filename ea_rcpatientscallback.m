function ea_rcpatientscallback(handles)

if get(handles.recentpts,'Value')==1
    return
end

load([ea_getearoot,'common',filesep,'ea_recentpatients.mat']);
if iscell(fullrpts)
    fullrpts=fullrpts(get(handles.recentpts,'Value')-1);
end

if strcmp('No recent patients found',fullrpts)
   return
end

ea_load_pts(handles,fullrpts);

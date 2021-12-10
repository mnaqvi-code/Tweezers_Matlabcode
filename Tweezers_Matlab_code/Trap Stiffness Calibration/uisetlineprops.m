function varargout = uisetlineprops(varargin)
% UISETLINEPROPS is a GUI tool for setting line and marker properties
%    Usage:
%      uisetlineprops(H) where H is the handle to a line object
%
%    Example:
%      s=peaks;
%      handle=plot(s(:,25));
%      props=uisetlineprops(handle)

%  MJ  - Mar 7, 2007

% Begin initialization code - DO NOT EDIT
if nargin==0
  feval('help',mfilename); return
end
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uisetlineprops_OpeningFcn, ...
                   'gui_OutputFcn',  @uisetlineprops_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before uisetlineprops is made visible.
function uisetlineprops_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to uisetlineprops (see VARARGIN)

h=[];
choices.linestyles=get(handles.popupmenu_linestyle,'string');
choices.markerstyles=get(handles.popupmenu_marker,'string');
setappdata(hObject,'choices',choices);
obj.PropName={'visible','color','linestyle','linewidth','marker','markersize'};
obj.PropVal=getprops(handles);

if nargin>3
  h=varargin{1};  h=h(ishandle(h));
  objtype=get(h,'type'); isline=false(size(h));
  ii=strmatch('line',objtype,'exact'); isline(ii)=true;
  if any(~isline)
    disp('Warning: Only "line" objects are supported');    h=h(isline);
  end
  if isempty(h)
    disp('Warning: Invalid object handle.');   
  end
end

obj.handle=h;

if ~isempty(h)
  obj.PropVal=get(h(1),obj.PropName)'; rgb=round(obj.PropVal{2}.*255);
  if ~strcmp(obj.PropVal{5},'none')  obj.PropVal{5}=obj.PropVal{5}(1); end
  
  setprops(handles,obj.PropVal);
end

%keyboard
setappdata(hObject,'object',obj);


% Choose default command line output for uisetlineprops
handles.output = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes uisetlineprops wait for user response (see UIRESUME)
setrgbcolor(handles);
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = uisetlineprops_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
try
  varargout{1} = handles.output;
catch
  varargout{1} = [];
end

closereq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_Red_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'val',str2double(get(hObject,'str')));
function edit_Green_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'val',str2double(get(hObject,'str')));
function edit_Blue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'val',str2double(get(hObject,'str')));
function popupmenu_linestyle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_linewidth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'val',str2double(get(hObject,'str')));
function popupmenu_marker_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_markersize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'val',str2double(get(hObject,'str')));
function slider_linewidth_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%set(hObject,'userdata',get(hObject,'val'));
function slider_markersize_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%set(hObject,'userdata',get(hObject,'val'));
%%%%%%%%%%%%%%%------------------------------------------------
%%%%%%%%%%%%%%%------------------------------------------------
%%%%%%%%%%%%%%%------------------------------------------------
%%%%%%%%%%%%%%%------------------------------------------------


function checkbox_visibility_Callback(hObject, eventdata, handles)
obj=getappdata(handles.figure1,'object');
set(obj.handle,'visible',ooflag(get(hObject,'value')));


function edit_Red_Callback(hObject, eventdata, handles)
obj=getappdata(handles.figure1,'object');
str=get(hObject,'string');  num=str2double(str);
num0=get(hObject,'value');  str0=num2str(num0);
if num>255 | num<0, num=nan; end
if num>0 & num<1, num=num.*255; end
num=round(num); str=num2str(num);
if ~isnan(num)   
  set(hObject,'value',num,'string',str);
  clr=[num get(handles.edit_Green,'val') get(handles.edit_Blue,'val')]./255;
  sett(obj.handle,'color',clr);
  setrgbcolor(handles);
else
  set(hObject,'string',str0);
end


function edit_Green_Callback(hObject, eventdata, handles)
obj=getappdata(handles.figure1,'object');
str=get(hObject,'string');  num=str2double(str);
num0=get(hObject,'value');  str0=num2str(num0);
if num>255 | num<0, num=nan; end
if num>0 & num<1, num=num.*255; end
num=round(num);  str=num2str(num);
if ~isnan(num)   
  set(hObject,'value',num,'string',str);
  clr=[get(handles.edit_Red,'val') num get(handles.edit_Blue,'val')]./255;
  sett(obj.handle,'color',clr);
  setrgbcolor(handles);
else
  set(hObject,'string',str0);
end


function edit_Blue_Callback(hObject, eventdata, handles)
obj=getappdata(handles.figure1,'object');
str=get(hObject,'string');  num=str2double(str);
num0=get(hObject,'value');  str0=num2str(num0);
if num>255 | num<0, num=nan; end
if num>0 & num<1, num=num.*255; end
num=round(num);  str=num2str(num);
if ~isnan(num)   
  set(hObject,'value',num,'string',str);
  clr=[get(handles.edit_Red,'val') get(handles.edit_Green,'val') num]./255;
  sett(obj.handle,'color',clr);
  setrgbcolor(handles);
else
  set(hObject,'string',str0);
end


function pushbutton_uisetcolor_Callback(hObject, eventdata, handles)
obj=getappdata(handles.figure1,'object');
if isempty(obj.handle)
  clr=uisetcolor('Set Line Color');
else
  clr=uisetcolor(obj.handle,'Set Line Color');
end
rgb=round(clr.*255);
set(handles.edit_Red,'str',num2str(rgb(1)),'val',rgb(1));
set(handles.edit_Green,'str',num2str(rgb(2)),'val',rgb(2));
set(handles.edit_Blue,'str',num2str(rgb(3)),'val',rgb(3));
setrgbcolor(handles);


function popupmenu_linestyle_Callback(hObject, eventdata, handles)
obj=getappdata(handles.figure1,'object');
str=get(hObject,'string'); val=get(hObject,'value');
sett(obj.handle,'linestyle',str{val});

function edit_linewidth_Callback(hObject, eventdata, handles)
obj=getappdata(handles.figure1,'object');
str=get(hObject,'string');  num=str2double(str);
num0=get(hObject,'value');  str0=num2str(num0);
if num<0, num=1; end
%if num>0 & num<1, num=num.*255; end
str=num2str(num);
if ~isnan(num)   
  set(hObject,'value',num,'string',str);
  sett(obj.handle,'linewidth',num);
else
  set(hObject,'string',str0);
end

function slider_linewidth_Callback(hObject, eventdata, handles)
val=get(hObject,'value'); set(hObject,'value',2);
num=get(handles.edit_linewidth,'val');
if val>2  % up
  num=floor(num*4/3)*3/4+3/4;
else
  num=max([0 ceil(num*4/3)*3/4-3/4]);
end
set(handles.edit_linewidth,'val',num,'string',num2str(num));
obj=getappdata(handles.figure1,'object');
sett(obj.handle,'linewidth',num);


function popupmenu_marker_Callback(hObject, eventdata, handles)
obj=getappdata(handles.figure1,'object');
str=get(hObject,'string'); val=get(hObject,'value');
sett(obj.handle,'marker',str{val});
if strmatch(str{val},{'.','x','+','*'})
  ec=get(obj.handle,'MarkerEdgeColor');
  ac=get(get(obj.handle,'parent'),'Color');
  if strcmp(ec,'none') | isequal(ec,ac)
    set(obj.handle,'MarkerEdgeColor',1-ac);
  end
end

function edit_markersize_Callback(hObject, eventdata, handles)
obj=getappdata(handles.figure1,'object');
str=get(hObject,'string');  num=str2double(str);
num0=get(hObject,'value');  str0=num2str(num0);
if num<0, num=6; end
%if num>0 & num<1, num=num.*255; end
str=num2str(num);
if ~isnan(num)   
  set(hObject,'value',num,'string',str);
  sett(obj.handle,'markersize',num);
else
  set(hObject,'string',str0);
end

function slider_markersize_Callback(hObject, eventdata, handles)
val=get(hObject,'value'); set(hObject,'value',2);
num=get(handles.edit_markersize,'val');
if val>2  % up
  num=floor(num)+1;
else
  num=max([0 ceil(num)-1]);
end
set(handles.edit_markersize,'val',num,'string',num2str(num));
obj=getappdata(handles.figure1,'object');
sett(obj.handle,'markersize',num);


function pushbutton_MarkerFaceColor_Callback(hObject, eventdata, handles)
obj=getappdata(handles.figure1,'object');
%clr=get(obj.handle,'color'); clr0='';
%try  clr0=get(obj.handle,'MarkerFaceColor'); end
%if isnumeric(clr0), clr=clr0; end
clr=uisetcolor('Set Line Color');
if isequal(clr,0)
  set(obj.handle,'MarkerFaceColor','none');
else
  set(obj.handle,'MarkerFaceColor',clr);
end


function pushbutton_MarkerEdgeColor_Callback(hObject, eventdata, handles)
obj=getappdata(handles.figure1,'object');
%clr=get(obj.handle,'color'); clr0='';
%try  clr0=get(obj.handle,'MarkerEdgeColor'); end
%if isnumeric(clr0), clr=clr0; end
clr=uisetcolor('Set Line Color');
if isequal(clr,0)
  set(obj.handle,'MarkerEdgeColor','none');
else
  set(obj.handle,'MarkerEdgeColor',clr);
end


function pushbutton_OK_Callback(hObject, eventdata, handles)
obj=getappdata(handles.figure1,'object');
obj.PropVal=getprops(handles);  % get new properties
handles.output=cell2struct(obj.PropVal(:),obj.PropName(:));
guidata(handles.figure1,handles);
uiresume(handles.figure1);


function pushbutton_Cancel_Callback(hObject, eventdata, handles)
obj=getappdata(handles.figure1,'object');
%% set old properties
set(obj.handle,obj.PropName',obj.PropVal');
handles.output=cell2struct(obj.PropVal(:),obj.PropName(:));
guidata(handles.figure1,handles);
uiresume(handles.figure1);

function pushbutton_reset_Callback(hObject, eventdata, handles)
obj=getappdata(handles.figure1,'object');
%% set old properties
set(obj.handle,obj.PropName',obj.PropVal');
setprops(handles,obj.PropVal);
setrgbcolor(handles);
%handles.output=cell2struct(obj.PropVal(:),obj.PropName(:));
%guidata(handles.figure1,handles);


%%%%%(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
%%%%%(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
function setprops(handles,PropVal)
choices=getappdata(handles.figure1,'choices');
rgb=round(PropVal{2}.*255);
set(handles.checkbox_visibility,'val',ooflag(PropVal{1}));
set(handles.edit_Red,'str',num2str(rgb(1)),'val',rgb(1));
set(handles.edit_Green,'str',num2str(rgb(2)),'val',rgb(2));
set(handles.edit_Blue,'str',num2str(rgb(3)),'val',rgb(3));
set(handles.popupmenu_linestyle,'value',...
  strmatch(PropVal{3},choices.linestyles,'exact'));
set(handles.edit_linewidth,'str',num2str(PropVal{4}),...
  'val',PropVal{4});
set(handles.popupmenu_marker,'value',...
  strmatch(PropVal{5},choices.markerstyles,'exact'));
set(handles.edit_markersize,'str',num2str(PropVal{6}),...
  'val',PropVal{6});

function PropVal=getprops(handles)
%obj=getappdata(handles.figure1,'object');
choices=getappdata(handles.figure1,'choices');

PropVal{1}=ooflag(get(handles.checkbox_visibility,'val'));
rgb(1,1)=get(handles.edit_Red,'val');
rgb(1,2)=get(handles.edit_Green,'val');
rgb(1,3)=get(handles.edit_Blue,'val');
PropVal{2}=rgb./255;
PropVal{3}=choices.linestyles{get(handles.popupmenu_linestyle,'value')};
PropVal{4}=get(handles.edit_linewidth,'val');
PropVal{5}=choices.markerstyles{get(handles.popupmenu_marker,'value')};
PropVal{6}=get(handles.edit_markersize,'val');
%setappdata(handles.figure1,'object',obj);

function setrgbcolor(h)
rgb(1,1)=get(h.edit_Red,'val');
rgb(1,2)=get(h.edit_Green,'val');
rgb(1,3)=get(h.edit_Blue,'val');
clr=rgb./255;
set([h.textR h.textG h.textB],'BackGroundColor',clr);
set([h.textR h.textG h.textB],'ForeGroundColor',1-clr);


function out=ooflag(in)  %% convert 'on'/'off' to 0/1 and vice versa
oo={'off','on'};
if isnumeric(in)
  if isscalar(in)
    out=oo{logical(in)+1}; return   % return string 'on' or 'off'
  else
    out=oo(logical(in)+1); return   % return cell array {'on','off',...}
  end
elseif ischar(in)
  out=logical(strmatch(in,oo,'exact')-1); return  % 0, 1, or []
elseif iscell(in)
  out=false(size(in));  ii=strmatch(oo{2},in,'exact');
  out(ii)=true;    % return logical array [0 0 1 0 1 1 1 0 ...]
end

function sett(h,varargin)
h=h(ishandle(h));
try, set(h,varargin{:}); end




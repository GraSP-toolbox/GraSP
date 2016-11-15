function varargout = GFT_gui(varargin)
% GFT_GUI MATLAB code for GFT_gui.fig
%      GFT_GUI, by itself, creates a new GFT_GUI or raises the existing
%      singleton*.
%
%      H = GFT_GUI returns the handle to a new GFT_GUI or the handle to
%      the existing singleton*.
%
%      GFT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GFT_GUI.M with the given input arguments.
%
%      GFT_GUI('Property','Value',...) creates a new GFT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GFT_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GFT_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GFT_gui

% Last Modified by GUIDE v2.5 21-Sep-2015 16:31:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GFT_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @GFT_gui_OutputFcn, ...
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


% --- Executes just before GFT_gui is made visible.
function GFT_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GFT_gui (see VARARGIN)
handles.G = varargin{1};
handles.f = varargin{2};
if numel(varargin) > 2
    handles.paramplot = varargin{3};
else
    handles.paramplot = struct();
end
handles.G = gsp_compute_fourier_basis(handles.G);
handles.G.f = handles.f;
handles.G.f_tilde = handles.G.U' * handles.f;

set(handles.figure1, 'Renderer', 'OpenGL');

% Flags
handles.eventoSetado = 1;
handles.redline = 0;
handles.x = -1;
handles.flag = 0;
handles.redline3 = 0;
handles.x3 = -1;
handles.flag3 = 0;

% Setting the signal plot
axes(handles.axes1);
handles.f_fig = stem(handles.f, 'Color',[0 0.447058826684952 0.74117648601532]);

% Setting the graph plot
% handles.paramplot.show_edges = 1;
% handles.paramplot.bar = 1;
% handles.paramplot.cp = [0 0 0];
% Do not return a handle
guidata(hObject, handles);
updateGraphPlot(hObject,handles);

%Setting the GFT plot
axes(handles.axes3);
handles.gft = handles.G.f_tilde;
handles.gft_fig = stem(handles.gft, 'Color',[0 0.447058826684952 0.74117648601532]);


% Setting the default callbacks
set(handles.figure1,'WindowButtonDownFcn',{@figure1_WindowButtonDownFcn,hObject});
set(handles.figure1,'WindowButtonUpFcn',{@figure1_WindowButtonUpFcn,hObject});
set(handles.figure1,'WindowButtonMotionFcn',{@figure1_WindowButtonMotionFcn,hObject});


%set(handles.axes1,'xlim',[-10 10],'ylim',[-10 10]);
%handles.redbox=patch([1, 2, 2, 1, 1],[2, 2, 3, 3, 2], 'r');
% Choose default command line output for GFT_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
updateIndex(hObject, handles);

% UIWAIT makes GFT_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GFT_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
set(handles.figure1,'WindowButtonMotionFcn',{@my_MouseMoveFcn,hObject});
guidata(hObject,handles)

function figure3_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
set(handles.figure1,'WindowButtonMotionFcn',{@my_MouseMoveFcn3,hObject});
guidata(hObject,handles)

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
set(handles.figure1,'WindowButtonMotionFcn',{@figure1_WindowButtonMotionFcn,hObject});

handles.gft = handles.G.U' * handles.f;
%updateGFTPlot(hObject,handles);
axes(handles.axes3);
delete(handles.gft_fig);
handles.gft_fig = stem(handles.gft, 'Color',[0 0.447058826684952 0.74117648601532]);


updateGraphPlot(hObject,handles);

guidata(hObject,handles);

function figure3_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
set(handles.figure1,'WindowButtonMotionFcn',{@figure3_WindowButtonMotionFcn,hObject});

% delete(handles.f_fig);
% updateSignalPlot(hObject,handles);
handles.f = handles.G.U * handles.gft;
axes(handles.axes1);
delete(handles.f_fig);
handles.f_fig = stem(handles.f, 'Color',[0 0.447058826684952 0.74117648601532]);

updateGraphPlot(hObject,handles);

guidata(hObject,handles);

function my_MouseMoveFcn(hObject, eventdata, handles)
handles=guidata(hObject);
axes(handles.axes1);
mousepos=get(handles.axes1,'CurrentPoint');
set(handles.text2, 'String', mousepos(1,1));
set(handles.text3, 'String', mousepos(1,2));
handles.f(handles.x) = mousepos(1,2);

if (handles.flag == 1)
    delete(handles.redline);
    handles.flag = 0;
end
delete(handles.f_fig);
handles.f_fig = stem(handles.f, 'Color',[0 0.447058826684952 0.74117648601532]);

handles.gft = handles.G.U' * handles.f;
%updateGFTPlot(hObject,handles);
axes(handles.axes3);
delete(handles.gft_fig);
handles.gft_fig = stem(handles.gft, 'Color',[0 0.447058826684952 0.74117648601532]);

updateGraphPlot(hObject,handles);

guidata(hObject, handles);
updateIndex(hObject, handles);
% x=get(handles.redbox,'xdata');
% y=get(handles.redbox,'ydata');
% x_rel=x-x(1);
% y_rel=y-y(1);
% set(handles.redbox,'xdata',mousepos(1,1)+x_rel);
% set(handles.redbox,'ydata',mousepos(1,2)+y_rel);

%function my_MouseMoveFcn3(obj,event,hObject)
function my_MouseMoveFcn3(hObject, eventdata, handles)
handles=guidata(hObject);
axes(handles.axes3);
mousepos=get(handles.axes3,'CurrentPoint');
set(handles.text2, 'String', mousepos(1,1));
set(handles.text3, 'String', mousepos(1,2));
handles.gft(handles.x3) = mousepos(1,2);

if (handles.flag3 == 1)
    delete(handles.redline3);
    handles.flag3 = 0;
end
delete(handles.gft_fig);
handles.gft_fig = stem(handles.gft, 'Color',[0 0.447058826684952 0.74117648601532]);

handles.f = handles.G.U * handles.gft;
%updateSignalPlot(hObject, handles);
axes(handles.axes1);
delete(handles.f_fig);
handles.f_fig = stem(handles.f, 'Color',[0 0.447058826684952 0.74117648601532]);

updateGraphPlot(hObject, handles);

guidata(hObject, handles);
updateIndex(hObject, handles);

% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
handles=guidata(hObject);
axes(handles.axes1);
hold on;
mousepos=get(handles.axes1,'CurrentPoint');
mousepos3 = get(handles.axes3,'CurrentPoint');
set(handles.text2, 'String', mousepos(1,1));
set(handles.text3, 'String', mousepos(1,2));

% If outside
if mousepos(1,2) > handles.axes1.YLim(2) || mousepos(1,2) < handles.axes1.YLim(1) || mousepos(1,1) > length(handles.f) || mousepos(1,1) < 0
    if mousepos3(1,1) > -0.1
        set(handles.figure1,'WindowButtonMotionFcn',{@figure3_WindowButtonMotionFcn,hObject});
        %set(handles.figure1,'WindowButtonMotionFcn',{@my_MouseMoveFcn3,hObject});
    else
        if handles.eventoSetado == 1
            set(handles.figure1,'WindowButtonDownFcn','');
            set(handles.figure1,'WindowButtonUpFcn','');
            handles.eventoSetado = 0;
        end
    end
    if (handles.flag == 1)
        handles.flag = 0;
        delete(handles.redline);
        guidata(hObject,handles);
    end 
    return
end

if (round(mousepos(1,1)) == handles.x) || (round(mousepos(1,1)) > length(handles.f)) || (round(mousepos(1,1)) < 1)
    return
end

if (handles.flag == 1)
    delete(handles.redline);
end
set(handles.figure1,'WindowButtonDownFcn',{@figure1_WindowButtonDownFcn,hObject});
set(handles.figure1,'WindowButtonUpFcn',{@figure1_WindowButtonUpFcn,hObject});
handles.eventoSetado = 1;
%set(handles.figure1,'WindowButtonMotionFcn',{@my_MouseMoveFcn,hObject});
handles.flag = 1;
handles.x = round(mousepos(1,1));
handles.redline = stem(handles.x,handles.f(handles.x), 'r');
guidata(hObject,handles);

% --- Executes on mouse motion over figure - except title and menu.
function figure3_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
handles=guidata(hObject);
axes(handles.axes3);
hold on;
mousepos=get(handles.axes3,'CurrentPoint');
set(handles.text2, 'String', mousepos(1,1));
set(handles.text3, 'String', mousepos(1,2));

if mousepos(1,2) > handles.axes3.YLim(2) || mousepos(1,2) < handles.axes3.YLim(1) || mousepos(1,1) > length(handles.gft) || mousepos(1,1) < 0
    if mousepos(1,1) < -0.1
        set(handles.figure1,'WindowButtonMotionFcn',{@figure1_WindowButtonMotionFcn,hObject});
        %set(handles.figure1,'WindowButtonMotionFcn',{@my_MouseMoveFcn,hObject});
    else
        set(handles.figure1,'WindowButtonDownFcn','');
        set(handles.figure1,'WindowButtonUpFcn','');
        %set(handles.figure1,'WindowButtonMotionFcn','');
    end
    if (handles.flag3 == 1)
        handles.flag3 = 0;
        delete(handles.redline3);
        guidata(hObject,handles);
    end 
    return
end

if (round(mousepos(1,1)) == handles.x3) || (round(mousepos(1,1)) > length(handles.gft)) || (round(mousepos(1,1)) < 1)
    return
end
if (handles.flag3 == 1)
    delete(handles.redline3);
end
set(handles.figure1,'WindowButtonDownFcn',{@figure3_WindowButtonDownFcn,hObject});
set(handles.figure1,'WindowButtonUpFcn',{@figure3_WindowButtonUpFcn,hObject});
%set(handles.figure1,'WindowButtonMotionFcn',{@my_MouseMoveFcn3,hObject});
handles.flag3 = 1;
handles.x3 = round(mousepos(1,1));
% if exist('handles.redline')
%     delete(handles.redline);
% end
handles.redline3 = stem(handles.x3,handles.gft(handles.x3), 'r');
guidata(hObject,handles);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles.axes1
handles.f

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
eigenVector = str2num(get(handles.edit3, 'String'));

if isempty(eigenVector)
    return
end
temp = zeros(size(handles.gft));
temp(eigenVector) = 1; %handles.gft(eigenVector);
handles.gft = temp;

handles.f = handles.G.U * handles.gft;
% Replotting all three graphs
%updateAllPlots(hObject,handles);
axes(handles.axes3);
delete(handles.gft_fig);
handles.gft_fig = stem(handles.gft, 'Color',[0 0.447058826684952 0.74117648601532]);

axes(handles.axes1);
delete(handles.f_fig);
handles.f_fig = stem(handles.f, 'Color',[0 0.447058826684952 0.74117648601532]);

axes(handles.axes2);
gsp_plot_signal(handles.G,handles.f, handles.paramplot);
%updateGraphPlot(hObject,handles);

guidata(hObject,handles);
updateIndex(hObject, handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
lowCutOff = str2num(get(handles.edit1, 'String'));
if isempty(lowCutOff)
    return
end
handles.gft(lowCutOff+1:end) = 0;

handles.f = handles.G.U * handles.gft;
% Replotting all three graphs
axes(handles.axes3);
delete(handles.gft_fig);
handles.gft_fig = stem(handles.gft, 'Color',[0 0.447058826684952 0.74117648601532]);

axes(handles.axes1);
delete(handles.f_fig);
handles.f_fig = stem(handles.f, 'Color',[0 0.447058826684952 0.74117648601532]);

axes(handles.axes2);
gsp_plot_signal(handles.G,handles.f, handles.paramplot);

%updateAllPlots(hObject,handles);
guidata(hObject,handles);
updateIndex(hObject, handles);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
highCutOff = str2num(get(handles.edit2, 'String'));
if isempty(highCutOff)
    return
end
handles.gft(1:highCutOff) = 0;

handles.f = handles.G.U * handles.gft;
% Replotting all three graphs
axes(handles.axes3);
delete(handles.gft_fig);
handles.gft_fig = stem(handles.gft, 'Color',[0 0.447058826684952 0.74117648601532]);

axes(handles.axes1);
delete(handles.f_fig);
handles.f_fig = stem(handles.f, 'Color',[0 0.447058826684952 0.74117648601532]);

axes(handles.axes2);
gsp_plot_signal(handles.G,handles.f, handles.paramplot);
%updateAllPlots(hObject,handles);
guidata(hObject,handles);
updateIndex(hObject, handles);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

% Recovering original coeficients;
handles.gft = handles.G.f_tilde;
% Recalculating the original signal;
handles.f = handles.G.U * handles.gft;
% Replotting all three graphs
axes(handles.axes3);
delete(handles.gft_fig);
handles.gft_fig = stem(handles.gft, 'Color',[0 0.447058826684952 0.74117648601532]);

axes(handles.axes1);
delete(handles.f_fig);
handles.f_fig = stem(handles.f, 'Color',[0 0.447058826684952 0.74117648601532]);

axes(handles.axes2);
gsp_plot_signal(handles.G,handles.f, handles.paramplot);

% handles.f = handles.G.U * handles.gft;
% axes(handles.axes1);
% delete(handles.f_fig);
% handles.f_fig = stem(handles.f, 'Color',[0 0.447058826684952 0.74117648601532]);
%updateAllPlots(hObject, handles)

guidata(hObject,handles);
updateIndex(hObject, handles);

function updateSignalPlot(hObject,handles)
handles=guidata(hObject);
handles.f = handles.G.U * handles.gft;
axes(handles.axes1);
delete(handles.f_fig);
handles.f_fig = stem(handles.f, 'Color',[0 0.447058826684952 0.74117648601532]);

guidata(hObject,handles);

function updateGraphPlot(hObject,handles)
handles=guidata(hObject);
axes(handles.axes2);
gsp_plot_signal(handles.G,handles.f, handles.paramplot);
%updateIndex(hObject, handles);
guidata(hObject,handles);

function updateGFTPlot(hObject,handles)
handles=guidata(hObject);
axes(handles.axes3);
delete(handles.gft_fig);
handles.gft_fig = stem(handles.gft, 'Color',[0 0.447058826684952 0.74117648601532]);

guidata(hObject,handles);

function updateAllPlots(hObject, handles)
updateSignalPlot(hObject,handles);
updateGraphPlot(hObject,handles);
updateGFTPlot(hObject,handles)

function updateIndex(hObject, handles)
handles=guidata(hObject);
[indiceAbsoluto, indiceEnergia] = graphFilter(handles.gft);
set(handles.index, 'String',indiceEnergia);
guidata(hObject,handles);

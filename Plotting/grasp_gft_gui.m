%MATLAB code for grasp_gft_gui.fig. Displays a GUI showing a graph signal
%in the direct (vertex) domain and the Fourier domain. Ideal filtering can
%be performed interactively.
%
%   GRASP_GFT_GUI(graph, signal) creates the GUI for the given graph and
%   signal.
%
%   GRASP_GFT_GUI(..., options) also provides options to display the graph
%   (see GRASP_SHOW_GRAPH for available options). options can either be
%   successive pair of arguments (name, value) or a matlab structure.
%
%   H = GRASP_GFT_GUI() returns the handle to a new GRASP_GFT_GUI or the
%   handle to the existing singleton.
%
%   GRASP_GFT_GUI('CALLBACK', hObject, eventData, handles,...) calls the
%   local function named CALLBACK in GRASP_GFT_GUI.m with the given input
%   arguments.
%
%   GRASP_GFT_GUI('Property','Value',...) creates a new GRASP_GFT_GUI or
%   raises the existing singleton.  Starting from the left, property value
%   pairs are applied to the GUI before grasp_gft_gui_OpeningFcn gets
%   called.  An unrecognized property name or invalid value makes property
%   application stop. All inputs are passed to grasp_gft_gui_OpeningFcn via
%   varargin.
%
% Authors:
%  - Diego R. C. Silva <diego at ect dot ufrn dot br>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Diego R. C. Silva, Universidade Federal do Rio Grande do Norte,
% BRAZIL (2015)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2016-2018)
% 
% diego at ect dot ufrn dot br
% benjamin.girault@usc.edu
% 
% This software is a computer program whose purpose is to provide a Matlab
% / Octave toolbox for handling and displaying graph signals.
% 
% This software is governed by the CeCILL license under French law and
% abiding by the rules of distribution of free software.  You can  use, 
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info". 
% 
% As a counterpart to the access to the source code and  rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty  and the software's author,  the holder of the
% economic rights,  and the successive licensors  have only  limited
% liability. 
% 
% In this respect, the user's attention is drawn to the risks associated
% with loading,  using,  modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean  that it is complicated to manipulate,  and  that  also
% therefore means  that it is reserved for developers  and  experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or 
% data to be ensured and,  more generally, to use and operate it in the 
% same conditions as regards security. 
% 
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL license and that you accept its terms.

function varargout = grasp_gft_gui(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @grasp_gft_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @grasp_gft_gui_OutputFcn, ...
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


% --- Executes just before grasp_gft_gui is made visible.
function grasp_gft_gui_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to grasp_gft_gui (see VARARGIN)

handles.graph  = varargin{1};
handles.orig_signal = varargin{2};
if numel(varargin) > 2
    if numel(varargin) > 3
        handles.paramplot = cell2struct(varargin(4:2:end), varargin(3:2:end), 2);
    else
        handles.paramplot = varargin{3:end};
    end
else
    handles.paramplot = struct();
end
if handles.graph.Finv == 0
    error('Please compute the graph Fourier transform first using grasp_eigendecomposition!');
end
handles.orig_signal_hat = grasp_fourier(handles.graph, handles.orig_signal);

% Some data required within handles
handles.redline1 = 0;
handles.vertex_id = -1;
handles.redlineflag1 = 0;
handles.redline3 = 0;
handles.freq_id = -1;
handles.redlineflag3 = 0;

% Some aliases
handles.graph_ax = handles.axes2;
handles.signal_ax = handles.axes1;
handles.gft_ax = handles.axes3;

% Setting the signal plot
handles.signal = handles.orig_signal;

cla(handles.signal_ax);
handles.signal_fig = stem(handles.signal_ax, handles.signal, 'Color', [0 0.447058826684952 0.74117648601532]);
xlim(handles.signal_ax, [0.5 grasp_nb_nodes(handles.graph) + 0.5]);
xlabel(handles.signal_ax, 'Vertex');
hold(handles.signal_ax, 'on');

% Setting the graph plot
guidata(hObject, handles);
cla(handles.graph_ax);
options.node_values = handles.signal;
options = grasp_merge_structs(handles.paramplot, options);
[handles.nodes_h, handles.edges_h] = grasp_show_graph(handles.graph_ax, handles.graph, options);

%Setting the GFT plot
handles.signal_hat = handles.orig_signal_hat;
handles.graph.freqs = grasp_frequencies(handles.graph);

cla(handles.gft_ax);
handles.gft_fig = stem(handles.gft_ax, handles.graph.freqs, handles.signal_hat, 'Color',[0 0.447058826684952 0.74117648601532]);
xlim(handles.gft_ax, [0 0.5]);
xlabel('Graph Frequency');
hold(handles.gft_ax, 'on');

% Setting the default callbacks
guidata(hObject, handles);
set(handles.figure1, 'WindowButtonDownFcn', @(hobj, evt) figure1_WindowButtonDownFcn(hobj, evt, guidata(hobj)));
set(handles.figure1, 'WindowButtonUpFcn', @(hobj, evt) figure1_WindowButtonUpFcn(hobj, evt, guidata(hobj)));
set(handles.figure1, 'WindowButtonMotionFcn', @(hobj, evt) figure1_WindowButtonMotionFcn(hobj, evt, guidata(hobj)));

% Choose default command line output for grasp_gft_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes grasp_gft_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = grasp_gft_gui_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, ~, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.redlineflag1
    set(handles.figure1, 'WindowButtonMotionFcn', @(hobj, evt) my_MouseMoveFcn1(hobj, evt, guidata(hobj)));
elseif handles.redlineflag3
    set(handles.figure1, 'WindowButtonMotionFcn', @(hobj, evt) my_MouseMoveFcn3(hobj, evt, guidata(hobj)));
end
guidata(hObject,handles)

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, ~, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1, 'WindowButtonMotionFcn', @(hobj, evt) figure1_WindowButtonMotionFcn(hobj, evt, guidata(hobj)));
guidata(hObject, handles);

function my_MouseMoveFcn1(hObject, ~, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mousepos = get(handles.signal_ax, 'CurrentPoint');
set(handles.text2, 'String', mousepos(1, 1));
set(handles.text3, 'String', mousepos(1, 2));
handles.signal(handles.vertex_id) = mousepos(1, 2);
handles.signal_hat = grasp_fourier(handles.graph, handles.signal);

if (handles.redlineflag1 == 1)
    delete(handles.redline1);
    handles.redlineflag1 = 0;
end

updateAllPlots(handles);

guidata(hObject, handles);

function my_MouseMoveFcn3(hObject, ~, handles)
% hObject    handle to figure3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mousepos = get(handles.gft_ax, 'CurrentPoint');
set(handles.text2, 'String', mousepos(1, 1));
set(handles.text3, 'String', mousepos(1, 2));
handles.signal_hat(handles.freq_id) = mousepos(1, 2);
handles.signal = grasp_fourier_inverse(handles.graph, handles.signal_hat);

if (handles.redlineflag3 == 1)
    delete(handles.redline3);
    handles.redlineflag3 = 0;
end

updateAllPlots(handles);

guidata(hObject, handles);

% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, ~, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mousepos_signal = get(handles.signal_ax, 'CurrentPoint');
mousepos_gft = get(handles.gft_ax, 'CurrentPoint');
set(handles.text2, 'String', mousepos_signal(1, 1));
set(handles.text3, 'String', mousepos_signal(1, 2));

% Where is the cursor?
mouse_in_signal_axes = mousepos_signal(1, 1) >= handles.signal_ax.XLim(1) && mousepos_signal(1, 1) <= handles.signal_ax.XLim(2);
mouse_in_signal_axes = mouse_in_signal_axes && mousepos_signal(1, 2) >= handles.signal_ax.YLim(1) && mousepos_signal(1, 2) <= handles.signal_ax.YLim(2);

mouse_in_gft_axes = mousepos_gft(1, 1) >= handles.gft_ax.XLim(1) && mousepos_gft(1, 1) <= handles.gft_ax.XLim(2);
mouse_in_gft_axes = mouse_in_gft_axes && mousepos_gft(1, 2) >= handles.gft_ax.YLim(1) && mousepos_gft(1, 2) <= handles.gft_ax.YLim(2);

% Are we outside some axes, and have previously set a red stem? => delete
% it
if ~mouse_in_signal_axes && handles.redlineflag1
    handles.redlineflag1 = 0;
    handles.vertex_id = -1;
    delete(handles.redline1);
    guidata(hObject, handles);
end

if ~mouse_in_gft_axes && handles.redlineflag3
    handles.redlineflag3 = 0;
    handles.freq_id = -1;
    delete(handles.redline3);
    guidata(hObject, handles);
end

% Outside of both axes? finish (nothing to update)
if ~mouse_in_signal_axes && ~mouse_in_gft_axes
    return
end

% Within gft axes bounds? get the closest stem
closest_freq_id = -1;
if mouse_in_gft_axes
    gft_stem_points = [handles.gft_fig.XData' handles.gft_fig.YData'];
    distances = sum((gft_stem_points - repmat(mousepos_gft(1, 1:2), size(gft_stem_points, 1), 1)) .^ 2, 2);
    [~, closest_freq_id] = min(distances);
end

% Have we set a redline and are still within its bounds?
if handles.redlineflag1
    within_bounds = round(mousepos_signal(1, 1)) == handles.vertex_id;
    if within_bounds
        return;
    end
end

if handles.redlineflag3
    if closest_freq_id == handles.freq_id
        return;
    end
end

% If a red stem remains, it means it changed => deleting before update
if handles.redlineflag1
    delete(handles.redline1);
end
if handles.redlineflag3
    delete(handles.redline3);
end

% New red stem
if mouse_in_signal_axes
    handles.redlineflag1 = 1;
    handles.vertex_id = round(mousepos_signal(1, 1));
    handles.redline1 = stem(handles.signal_ax, handles.vertex_id, handles.signal(handles.vertex_id), 'r');
end
if mouse_in_gft_axes
    handles.redlineflag3 = 1;
    handles.freq_id = closest_freq_id;
    handles.redline3 = stem(handles.gft_ax, handles.graph.freqs(closest_freq_id), handles.signal_hat(closest_freq_id), 'r');
end

guidata(hObject, handles);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp(handles.signal)

function edit1_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end


function edit2_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end



function edit3_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

eigenVector = floor(str2double(get(handles.edit3, 'String')));

if isnan(eigenVector)
    return
end
handles.signal_hat = zeros(size(handles.signal_hat));
handles.signal_hat(eigenVector) = 1;

handles.signal = grasp_fourier_inverse(handles.graph, handles.signal_hat);

% Replotting all three graphs
updateAllPlots(handles);

guidata(hObject, handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles=guidata(hObject);

lowCutOff = str2double(get(handles.edit1, 'String'));
if isnan(lowCutOff)
    return
end
handles.signal_hat(grasp_frequencies(handles.graph) > lowCutOff) = 0;
handles.signal = grasp_fourier_inverse(handles.graph, handles.signal_hat);

% Replotting all three graphs
updateAllPlots(handles);

guidata(hObject, handles);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

highCutOff = str2double(get(handles.edit2, 'String'));
if isnan(highCutOff)
    return
end
handles.signal_hat(grasp_frequencies(handles.graph) < highCutOff) = 0;
handles.signal = grasp_fourier_inverse(handles.graph, handles.signal_hat);

% Replotting all three graphs
updateAllPlots(handles);

guidata(hObject, handles);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Recovering original coeficients;
handles.signal_hat = handles.orig_signal_hat;
handles.signal = handles.orig_signal;

updateAllPlots(handles);
guidata(hObject, handles);


% -----------------------------------------------------
% Helper functions

function updateSignalPlot(handles)
set(handles.signal_fig, 'YData', handles.signal);
refreshdata(handles.signal_fig);

function updateGraphPlot(handles)
grasp_scatter_update(handles.nodes_h, handles.signal);
refreshdata(handles.nodes_h);

function updateGFTPlot(handles)
set(handles.gft_fig, 'YData', handles.signal_hat);
refreshdata(handles.gft_fig);

function updateAllPlots(handles)
updateSignalPlot(handles);
updateGraphPlot(handles);
updateGFTPlot(handles);

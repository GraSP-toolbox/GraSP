%MATLAB code for grasp_animate_iterated_operator_gui.fig. Iterates an 
%operator on a signal and shows an animation of its output.
%
%   GRASP_ANIMATE_ITERATED_OPERATOR_GUI(graph, operator, signal) creates
%   the GUI for the given graph and signal. The operator is then iterated
%   on the signal.
%
%   GRASP_ANIMATE_ITERATED_OPERATOR_GUI(..., show_real_imag) boolean on
%   whether to show the real and imaginary part of the output signal on
%   separate axes.
%
%   GRASP_ANIMATE_ITERATED_OPERATOR_GUI(..., file_basename) base filename
%   for the snapshots.
%
%   H = GRASP_ANIMATE_ITERATED_OPERATOR_GUI() returns the handle to a new
%   GRASP_ANIMATE_ITERATED_OPERATOR_GUI or the handle to the existing
%   singleton.
%
%   GRASP_ANIMATE_ITERATED_OPERATOR_GUI('CALLBACK', hObject, eventData, handles,...)
%   calls the local function named CALLBACK in
%   GRASP_ANIMATE_ITERATED_OPERATOR_GUI.m with the given input arguments.
%
%   GRASP_ANIMATE_ITERATED_OPERATOR_GUI('Property','Value',...) creates a
%   new GRASP_ANIMATE_ITERATED_OPERATOR_GUI or raises the existing
%   singleton.  Starting from the left, property value pairs are applied to
%   the GUI before grasp_animate_iterated_operator_gui_OpeningFcn gets
%   called.  An unrecognized property name or invalid value makes property
%   application stop. All inputs are passed to
%   grasp_animate_iterated_operator_gui_OpeningFcn via varargin.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2016)
% 
% benjamin.girault@ens-lyon.fr
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

function varargout = grasp_animate_iterated_operator_gui(varargin)


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @grasp_animate_iterated_operator_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @grasp_animate_iterated_operator_gui_OutputFcn, ...
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



function [data, s, sraw] = apply_operator(data)
    data.current_iteration = round(data.current_iteration * data.nb_non_int) / data.nb_non_int;
    if data.prev_int_power_ind ~= floor(data.current_iteration)
        data.prev_int_power_ind = floor(data.current_iteration);
        data.prev_int_power = data.operator ^ data.prev_int_power_ind;
    end
    cur_iteration_non_int_ind = round((data.current_iteration - data.prev_int_power_ind) / abs(data.step) + 1);
    if cur_iteration_non_int_ind == 1
        sraw = data.prev_int_power * data.signal;
    else
        sraw = (data.prev_int_power * data.iteration_non_int{cur_iteration_non_int_ind}) * data.signal;
    end

    switch data.value
        case 1
            s = abs(sraw);
        case 2
            s = abs(sraw) .^ 2;
        case 3
            s = angle(sraw);
        case 4
            s = real(sraw);
        case 5
            s = abs(real(sraw));
        case 6
            s = imag(sraw);
        case 7
            s = abs(imag(sraw));
    end
% end

function main_loop_body(obj, ~)
    timer_data = get(obj, 'UserData');
    
    timer_data.current_iteration = timer_data.current_iteration + timer_data.step;
    [timer_data, s, sraw] = apply_operator(timer_data);
    grasp_scatter_update(timer_data.nodes_handle, s);
    refreshdata(timer_data.nodes_handle);
    if timer_data.show_real_imag ~= 0
        grasp_scatter_update(timer_data.nodes_handle_real,real(sraw));
        grasp_scatter_update(timer_data.nodes_handle_imag, imag(sraw));
        refreshdata(timer_data.nodes_handle_real);
        refreshdata(timer_data.nodes_handle_imag);
    end
    set(timer_data.iteration_text_handle, 'String', sprintf('%c=%f', char(964), timer_data.current_iteration));
    set(obj, 'UserData', timer_data);
    drawnow;
% end

function precompute(obj, ~)
    % Data
    timer_data = get(obj, 'UserData');
    
    % Precomputations
    for i = 1:timer_data.nb_non_int
        timer_data.iteration_non_int{i} = timer_data.operator ^ (abs(timer_data.step) * (i - 1));
    end
    
    % Store data
    set(obj, 'UserData', timer_data);

    % Enable reset
    set(timer_data.reset_handle, 'Enable', 'on');

    % Enable animation
    set(timer_data.start_handle, 'Enable', 'on');
    set(timer_data.step_handle, 'Enable', 'on');
%end


% --- Executes just before grasp_animate_iterated_operator_gui is made visible.
function grasp_animate_iterated_operator_gui_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to grasp_animate_iterated_operator_gui (see VARARGIN)

% Choose default command line output for grasp_animate_iterated_operator_gui
handles.output = hObject;

% Input
file_basename = '/tmp/iterations';
show_real_imag = 0;
if numel(varargin) >= 4
    for i = 4:numel(varargin)
        if isnumeric(varargin{i})
            show_real_imag = varargin{i};
        else
            file_basename = varargin{i};
        end
    end
end
handles.file_basename = file_basename;
handles.timer_data = struct('graph', varargin{1},...
                            'operator', varargin{2},...
                            'signal', varargin{3} / norm(varargin{3}),...
                            'show_real_imag', show_real_imag,...
                            'start_handle', handles.pushbutton1,...
                            'reset_handle', handles.pushbutton3,...
                            'step_handle', handles.pushbutton5,...
                            'nodes_handle', 0,...
                            'nodes_handle_real', 0,...
                            'nodes_handle_imag', 0,...
                            'iteration_text_handle', handles.text4,...
                            'current_iteration', 0,...
                            'step', 0);

% Timer to start / pause / stop the animation
handles.timer = timer('ExecutionMode', 'fixedSpacing',...
                      'Period', 0.001,...
                      'TimerFcn', @main_loop_body,...
                      'UserData', handles.timer_data);
                  
% Timer for precomputations
handles.timer_precomp = timer('TimerFcn', @precompute, 'StartDelay', 0.1);

% UIWAIT makes grasp_animate_iterated_operator_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Axes
% if strcmp(displayed_values, 'angle')
%     colormap(handles.axes1, 'hsv');
% end
cla(handles.axes1);
hold(handles.axes1, 'on');
cla(handles.axes2);
hold(handles.axes2, 'on');
cla(handles.axes3);
hold(handles.axes3, 'on');
cmin = get(handles.slider1, 'Value');
cmax = get(handles.slider2, 'Value');
handles.timer_data.nodes_handle = grasp_show_graph(handles.axes1, handles.timer_data.graph, struct(...
                                             'value_scale', [cmin cmax],...
                                             'node_display_size', 1500));
if show_real_imag ~= 0
    handles.timer_data.nodes_handle_real = grasp_show_graph(handles.axes2, handles.timer_data.graph, struct(...
                                                      'value_scale', 0.2 * [-1 1],...
                                                      'node_display_size', 1000));
    grasp_set_blue_red_colormap(handles.axes2);
    handles.timer_data.nodes_handle_imag = grasp_show_graph(handles.axes3, handles.timer_data.graph, struct(...
                                                      'value_scale', 0.2 * [-1 1],...
                                                      'node_display_size', 1000));
    grasp_set_blue_red_colormap(handles.axes3);
end
switch get(handles.popupmenu1, 'Value')
    case 3
        colormap(handles.axes1, 'hsv');
    case 4
        grasp_set_blue_red_colormap(handles.axes1);
    case 6
        grasp_set_blue_red_colormap(handles.axes1);
    otherwise
        colormap(handles.axes1, 'jet');
end

set(handles.text4, 'String', sprintf('%c=%f', char(964), handles.timer_data.current_iteration));

% Enable / Disable buttons
set(handles.pushbutton1, 'Enable', 'off');      % Start
set(handles.pushbutton5, 'Enable', 'off');      % Step
set(handles.pushbutton3, 'Enable', 'on');       % Reset
set(handles.pushbutton2, 'Enable', 'off');      % Stop
set(handles.pushbutton6, 'Enable', 'on');       % Snapshot
set(handles.togglebutton1, 'Enable', 'off');    % Pause

% Enable parameters
set(handles.popupmenu1, 'Enable', 'on');
set(handles.edit1, 'Enable', 'on');
set(handles.slider1, 'Enable', 'on');
set(handles.slider2, 'Enable', 'on');

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = grasp_animate_iterated_operator_gui_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

% Force reset from user...
set(handles.pushbutton1, 'Enable', 'off');
set(handles.pushbutton5, 'Enable', 'off');

% Check sliders
value = get(hObject, 'Value');
if value ~= 3 || value ~=4 || value ~=6
    cmin = get(handles.slider1, 'Value');
    cmax = get(handles.slider2, 'Value');
    if cmin < 0
        cmin = 0;
    end
    if cmax < 0.1
        cmax = 0.1;
    end
    set(handles.slider1, 'Value', cmin);
    set(handles.slider2, 'Value', cmax);
    set(handles.text7, 'String', num2str(cmin));
    set(handles.text8, 'String', num2str(cmax));
    caxis(handles.axes1, [cmin cmax]);
end

% Angle is a special case
if value == 3
    colormap(handles.axes1, 'hsv');
    set(handles.slider1, 'Value', -1);
    set(handles.slider2, 'Value', 1);
    cmin = -1;
    cmax = 1;
    set(handles.text7, 'String', num2str(pi * cmin));
    set(handles.text8, 'String', num2str(pi * cmax));
    caxis(handles.axes1, pi * [cmin, cmax]);
elseif value == 4 || value == 6
    grasp_set_blue_red_colormap(handles.axes1);
else
    colormap(handles.axes1, 'default');
end


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1 (Start).
function pushbutton1_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Enable / Disable buttons
set(handles.pushbutton1, 'Enable', 'off');  % Start
set(handles.pushbutton5, 'Enable', 'off');  % Step
set(handles.pushbutton3, 'Enable', 'off');  % Reset
set(handles.pushbutton2, 'Enable', 'on');   % Stop
set(handles.pushbutton6, 'Enable', 'off');   % Snapshot
set(handles.togglebutton1, 'Enable', 'on'); % Pause
set(handles.togglebutton1, 'Value', get(handles.togglebutton1, 'Min'));

% Disable parameters
set(handles.popupmenu1, 'Enable', 'off');
set(handles.edit1, 'Enable', 'off');
set(handles.slider1, 'Enable', 'off');
set(handles.slider2, 'Enable', 'off');

% Get precomputations
timer_precomp_data = get(handles.timer_precomp, 'UserData');
timer_data = get(handles.timer, 'UserData');
timer_data.iteration_non_int = timer_precomp_data.iteration_non_int;

% Start the timer
set(handles.timer, 'ExecutionMode', 'fixedSpacing');
set(handles.timer, 'UserData', timer_data);
start(handles.timer);



% --- Executes on button press in pushbutton2 (Stop).
function pushbutton2_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Stop the timer
stop(handles.timer);

% Enable / Disable buttons
set(handles.pushbutton1, 'Enable', 'on');    % Start
set(handles.pushbutton5, 'Enable', 'on');    % Step
set(handles.pushbutton3, 'Enable', 'on');    % Reset
set(handles.pushbutton2, 'Enable', 'off');   % Stop
set(handles.pushbutton6, 'Enable', 'on');    % Snapshot
set(handles.togglebutton1, 'Enable', 'off'); % Pause

% Enable parameters
set(handles.popupmenu1, 'Enable', 'on');
set(handles.edit1, 'Enable', 'on');
set(handles.slider1, 'Enable', 'on');
set(handles.slider2, 'Enable', 'on');

% Reset the signal next step
handles.timer_data = get(handles.timer, 'UserData');
handles.timer_data.current_iteration = 0;
handles.timer_data.prev_int_power_ind = 0;
handles.timer_data.prev_int_power = eye(grasp_nb_nodes(handles.timer_data.graph));
set(handles.timer, 'UserData', handles.timer_data);

% Update figure
[handles.timer_data, s, sraw] = apply_operator(handles.timer_data);
% set(handles.timer_data.nodes_handle, 'CData', s);
grasp_scatter_update(handles.timer_data.nodes_handle, s);
refreshdata(handles.timer_data.nodes_handle);
if handles.timer_data.show_real_imag ~= 0
%     set(handles.timer_data.nodes_handle_real, 'CData', real(sraw));
    grasp_scatter_update(handles.timer_data.nodes_handle_real, real(sraw));
%     set(handles.timer_data.nodes_handle_imag, 'CData', imag(sraw));
    grasp_scatter_update(handles.timer_data.nodes_handle_imag, imag(sraw));
    refreshdata(handles.timer_data.nodes_handle_real);
    refreshdata(handles.timer_data.nodes_handle_imag);
end
set(handles.text4, 'String', sprintf('%c=%f', char(964), handles.timer_data.current_iteration));

% Update handles structure
guidata(hObject, handles);


function edit1_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% Force reset from user...
set(handles.pushbutton1, 'Enable', 'off');
set(handles.pushbutton5, 'Enable', 'off');



% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton1 (Pause).
function togglebutton1_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1

button_state = get(hObject,'Value');
if button_state == get(hObject,'Min')
	% Unpause
    set(handles.timer, 'ExecutionMode', 'fixedSpacing');
    start(handles.timer);
    
    % Disable / enable buttons
    set(handles.pushbutton5, 'Enable', 'off');  % Step
    set(handles.pushbutton6, 'Enable', 'off');  % Snapshot
    
    % Disable parameters
    set(handles.slider1, 'Enable', 'off');
    set(handles.slider2, 'Enable', 'off');
elseif button_state == get(hObject,'Max')
    % Pause
    stop(handles.timer);
    
    % Disable / enable buttons
    set(handles.pushbutton5, 'Enable', 'on');   % Step
    set(handles.pushbutton6, 'Enable', 'on');   % Snapshot
    
    % Enable parameters
    set(handles.slider1, 'Enable', 'on');
    set(handles.slider2, 'Enable', 'on');
end


% --- Executes on button press in pushbutton3 (Reset).
function pushbutton3_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Initializations
prev_step = handles.timer_data.step;
handles.timer_data.step = str2double(get(handles.edit1, 'String'));
if ~isfinite(handles.timer_data.step) || handles.timer_data.step == 0
    handles.timer_data.step = 0.1;
    set(handles.edit1, 'String', '0.1');
end
handles.timer_data.value = get(handles.popupmenu1, 'Value');

% Start precomputations
if prev_step ~= handles.timer_data.step
    % Disable Reset
    set(hObject, 'Enable', 'off');
    
    % Setup precomputations
    handles.timer_data.nb_non_int = abs(1 / handles.timer_data.step);
    if handles.timer_data.nb_non_int < 1
        handles.timer_data.nb_non_int = 1;
    end
    handles.timer_data.iteration_non_int = cell(handles.timer_data.nb_non_int, 1);
    
    % Start precomputations in background
    set(handles.timer_precomp, 'UserData', handles.timer_data);
    start(handles.timer_precomp);
else
    % Enable animation
    set(handles.timer_data.start_handle, 'Enable', 'on');
    set(handles.timer_data.step_handle, 'Enable', 'on');
end

% Current iteration
handles.timer_data.current_iteration = 0;
handles.timer_data.prev_int_power_ind = 0;
handles.timer_data.prev_int_power = eye(grasp_nb_nodes(handles.timer_data.graph));

% Sent data to timer
set(handles.timer, 'UserData', handles.timer_data);

% Update figure
[handles.timer_data, s, sraw] = apply_operator(handles.timer_data);
% set(handles.timer_data.nodes_handle, 'CData', s);
grasp_scatter_update(handles.timer_data.nodes_handle, s);
refreshdata(handles.timer_data.nodes_handle);
if handles.timer_data.show_real_imag ~= 0
%     set(handles.timer_data.nodes_handle_real, 'CData', real(sraw));
    grasp_scatter_update(handles.timer_data.nodes_handle_real, real(sraw));
%     set(handles.timer_data.nodes_handle_imag, 'CData', imag(sraw));
    grasp_scatter_update(handles.timer_data.nodes_handle_imag, imag(sraw));
    refreshdata(handles.timer_data.nodes_handle_real);
    refreshdata(handles.timer_data.nodes_handle_imag);
end
set(handles.text4, 'String', sprintf('%c=%f', char(964), handles.timer_data.current_iteration));

% Update handles structure
guidata(hObject, handles);


% --- Executes on slider movement.
function slider1_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

cmin = round(100 * get(hObject, 'Value')) / 100;
cmax = get(handles.slider2, 'Value');

if cmin > 0.9
    cmin = 0.9;
end
if cmin > cmax
    cmin = cmax - 0.1;
end
set(hObject, 'Value', cmin);
if get(handles.popupmenu1, 'Value') == 3
    set(handles.text7, 'String', num2str(pi * cmin));
    caxis(handles.axes1, pi * [cmin cmax]);
else
    set(handles.text7, 'String', num2str(cmin));
    caxis(handles.axes1, [cmin cmax]);
end


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

cmin = get(handles.slider1, 'Value');
cmax = round(100 * get(hObject, 'Value')) / 100;

if cmax < -0.9
    cmax = -0.9;
end
if cmin > cmax
    cmax = cmin + 0.1;
end
set(hObject, 'Value', cmax);
if get(handles.popupmenu1, 'Value') == 3
    set(handles.text8, 'String', num2str(pi * cmax));
    caxis(handles.axes1, pi * [cmin cmax]);
else
    set(handles.text8, 'String', num2str(cmax));
    caxis(handles.axes1, [cmin cmax]);
end


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton5 (Step).
function pushbutton5_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Wait for any non finished step job
if strcmp(get(handles.timer, 'ExecutionMode'), 'singleShot')
    wait(handles.timer);
end

% Get timer data
timer_data = get(handles.timer, 'UserData');

% Get precomputations
timer_precomp_data = get(handles.timer_precomp, 'UserData');
timer_data.iteration_non_int = timer_precomp_data.iteration_non_int;

% Set the timer
set(handles.timer, 'ExecutionMode', 'singleShot');
set(handles.timer, 'UserData', timer_data);

% Disable button
% set(handles.timer_data.step_handle, 'Enable', 'off');

% Update handles structure
guidata(hObject, handles);

% Start the timer
start(handles.timer);


% --- Executes on button press in pushbutton6 (Snapshot).
function pushbutton6_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

timer_data = get(handles.timer, 'UserData');
filename = sprintf('%s-%f.png', handles.file_basename, timer_data.current_iteration);
export_fig(filename, '-transparent', handles.axes1);

function varargout = prova(varargin)
% PROVA M-file for prova.fig
%      PROVA, by itself, creates a new PROVA or raises the existing
%      singleton*.
%
%      H = PROVA returns the handle to a new PROVA or the handle to
%      the existing singleton*.
%
%      PROVA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROVA.M with the given input arguments.
%
%      PROVA('Property','Value',...) creates a new PROVA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before prova_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to prova_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help prova

% Last Modified by GUIDE v2.5 02-Feb-2010 16:05:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prova_OpeningFcn, ...
                   'gui_OutputFcn',  @prova_OutputFcn, ...
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


% --- Executes just before prova is made visible.
function prova_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prova (see VARARGIN)

% Choose default command line output for prova
handles.output = hObject;

set(hObject,'toolbar','figure');
% Update handles structure
guidata(hObject, handles);


a2 = get(handles.a2,'Value');
a3 = get(handles.a3,'Value');
a1 = str2double(get(handles.a1,'String'));
a0 = str2double(get(handles.a0,'String'));
set(handles.txt_a3,'String',a3)
set(handles.txt_a2,'String',a2)
draw_shape(handles,a0,a1,a2,a3);

% UIWAIT makes prova wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prova_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function a0_Callback(hObject, eventdata, handles)
% hObject    handle to a0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a0 as text
%        str2double(get(hObject,'String')) returns contents of a0 as a double


% --- Executes during object creation, after setting all properties.
function a0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a1_Callback(hObject, eventdata, handles)
% hObject    handle to a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a1 as text
%        str2double(get(hObject,'String')) returns contents of a1 as a double


% --- Executes during object creation, after setting all properties.
function a1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function a2_Callback(hObject, eventdata, handles)
% hObject    handle to a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


a2 = get(hObject,'Value');
a3 = get(handles.a3,'Value');
a1 = str2double(get(handles.a1,'String'));
a0 = str2double(get(handles.a0,'String'));
set(handles.txt_a2,'String',a2)
draw_shape(handles,a0,a1,a2,a3);


% --- Executes during object creation, after setting all properties.
function a2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function a3_Callback(hObject, eventdata, handles)
% hObject    handle to a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
a2 = get(handles.a2,'Value');
a3 = get(hObject,'Value');
a1 = str2double(get(handles.a1,'String'));
a0 = str2double(get(handles.a0,'String'));
set(handles.txt_a3,'String',a3)
draw_shape(handles,a0,a1,a2,a3);

% --- Executes during object creation, after setting all properties.
function a3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over txt_a3.
function txt_a3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to txt_a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a2 = get(handles.a2,'Value');
a3 = get(handles.a3,'Value');
a1 = str2double(get(handles.a1,'String'));
a0 = str2double(get(handles.a0,'String'));
set(handles.txt_a3,'String',a3)
set(handles.txt_a2,'String',a2)
draw_shape(handles,a0,a1,a2,a3);


% --- Constraints and draw shape ---
function draw_shape(handles,a0,a1,a2,a3)

s= -1:0.001:1;
A= [ 1 0 0 0; 0 5 7 9; 0 20 42 72;0 60 210 504];
b = [-a0-a2;-a1-3*a3;-6*a3;-6*a3];

res =A\b;
a4 = res(1);
a5 = res(2);
a7 = res(3);
a9 = res(4);
    
x = a0 + a2*s.^2 + a4*s.^4;
y = a1.*s + a3*s.^3 + a5*s.^5 + a7*s.^7 + a9*s.^9;
x1 = -(a0 + a2*s.^2 + a4*s.^4);
z = zeros(1,length(x))
axes(handles.axes1);
cla;

plot(x,y)
hold on
plot(x1,y)
hold off
axis([-20 20 -20 20])

max_radius(handles,x,y);

save('shape', 'x', 'y','z')
%--------------------------------------

% --- Maximum radius ---
function max_radius(handles,x,y)
offset = 50;

P1=[x(1)^2 x(1) 1; x(1+offset)^2 x(1+offset) 1; x(1+2*offset)^2 x(1+2*offset) 1];
bp1 = [y(1); y(1+offset); y(1+2*offset)];
coef1 = P1\bp1;
ap1 = coef1(1);
bp1 = coef1(2);
cp1 = coef1(3);
if ap1 < 1e-10
    radius1 = inf;
else
    radius1 = 1/(2*ap1);
end

i = find(y==0);
P2=[y(i)^2 y(i) 1; y(i+offset)^2 y(i+offset) 1; y(i+2*offset)^2 y(i+2*offset) 1];
bp2 = [x(i); x(i+offset); x(i+2*offset)];
coef2 = P2\bp2;
ap2 = coef2(1);
bp2 = coef2(2);
cp2 = coef2(3);
if ap2 < 1e-10
    radius2 = inf;
else
    radius2 = 1/(2*ap2);
end

radius = min(radius1,radius2);
set(handles.max_rad,'String',radius)
%----------------------- 


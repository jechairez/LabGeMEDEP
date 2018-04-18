function varargout = Controllability_gui(varargin)
% CONTROLLABILITY_GUI MATLAB code for Controllability_gui.fig
%      CONTROLLABILITY_GUI, by itself, creates a new CONTROLLABILITY_GUI or raises the existing
%      singleton*.
%
%      H = CONTROLLABILITY_GUI returns the handle to a new CONTROLLABILITY_GUI or the handle to
%      the existing singleton*.
%
%      CONTROLLABILITY_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONTROLLABILITY_GUI.M with the given input arguments.
%
%      CONTROLLABILITY_GUI('Property','Value',...) creates a new CONTROLLABILITY_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Controllability_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Controllability_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Controllability_gui

% Last Modified by GUIDE v2.5 13-Dec-2017 16:48:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Controllability_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @Controllability_gui_OutputFcn, ...
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


% --- Executes just before Controllability_gui is made visible.
function Controllability_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Controllability_gui (see VARARGIN)

% Choose default command line output for Controllability_gui
handles.output = hObject;
global k ME MI MD MN MR MC MU MX nodes
k=2;
ME = lme(k); % equivalence matrix
MI = lmi(k); % inequivalence matrix 
MD = lmd(k); % OR matrix
MN = lmn(k); % NOT matrix
MR = lmr(k); % reduction power matrix
MC = lmc(k); % AND matrix
MU = lmu(k); % dummy matrix
MX = lm([2 1 1 2], 2); % xor
nodes={};
guidata(hObject, handles);

% UIWAIT makes Controllability_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Controllability_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb_browser.
function pb_browser_Callback(hObject, eventdata, handles)
% hObject    handle to pb_browser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fullpathname text A
[filename pathname]=uigetfile({'*.txt'},'File Selector');
fullpathname=strcat(pathname, filename);
text=fileread(fullpathname);
set(handles.txt_path2,'String',fullpathname);
set(handles.txt_content,'String',text);
set(handles.txt_content, 'HorizontalAlignment','left');
A=regexp(fileread(fullpathname),'\n','split');% Open archive for reading
guidata(hObject,handles);

% --- Executes on button press in pb_attractors.
function pb_attractors_Callback(hObject, eventdata, handles)
% hObject    handle to pb_attractors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global k nodes ME MI MD MN MR MC MU MX A
global Att_Land Att_Land_size
options = lmset('vars',nodes'); % Defining the vector of logical variables
%% Attractors and transition matrix of non-controlled network
[expr,vars]=GetAttractors(A,options,k,ME,MI,MD,MN,MR,MC,MU,MX); % Get attractors from uploaded txt
L=ctimes(expr{:}); % Get transition matrix from uploaded matrix from uploaded txt
[Att_Land,Att_Land_size]=ShowAttractors(L,k); % Get and show attractors //Check ShowAttractors function
set(handles.table_attractors,'Data',Att_Land);
set(handles.table_attractors,'ColumnName',nodes);
guidata(hObject,handles);



function et_nodes_Callback(hObject, eventdata, handles)
% hObject    handle to et_nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_nodes as text
%        str2double(get(hObject,'String')) returns contents of et_nodes as a double


% --- Executes during object creation, after setting all properties.
function et_nodes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_nodes.
function pb_nodes_Callback(hObject, eventdata, handles)
% hObject    handle to pb_nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global nodes
newvalue=get(handles.et_nodes,'String');
if isempty(nodes)
    nodes=newvalue;
    set(handles.txt_nodes,'String',nodes);
else
    nodes=vertcat(nodes,{newvalue});
    set(handles.txt_nodes,'String',nodes);
end
guidata(hObject,handles);



function et_steps_Callback(hObject, eventdata, handles)
% hObject    handle to et_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_steps as text
%        str2double(get(hObject,'String')) returns contents of et_steps as a double


% --- Executes during object creation, after setting all properties.
function et_steps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lb_connector.
function lb_connector_Callback(hObject, eventdata, handles)
% hObject    handle to lb_connector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb_connector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_connector


% --- Executes during object creation, after setting all properties.
function lb_connector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_connector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_controllability.
function pb_controllability_Callback(hObject, eventdata, handles)
global k nodes ME MI MD MN MR MC MU MX A
global Att_Land Att_Land_size n
global s m zz connector fullpathname
global Trajectories_a Trajectories_b U
m=1;
zz=1;
n=length(nodes)
s=str2double(get(handles.et_steps,'String'))
connector=get(handles.lb_connector,'Value');
if connector==1
    connectstr='MC'
else
    connectstr='MD'
end
options = lmset('vars',horzcat({'u1'},nodes')); % Defining the vector logical variables
for i=1:length(A) % Loop 1
    A{i}=strcat(A{i},[blanks(1) 'u1']); % Add control input on first gene , and then the next one and so on
    A{i}=sprintf('%s',[connectstr blanks(1) A{i}]); % Add logic operator on first gene, and then the next one and so on // can be OR, AND, etc operators
    [expr,vars]=GetAttractors(A,options,k,ME,MI,MD,MN,MR,MC,MU,MX); % Get attractors from modified txt (controlled network)
    L = ctimes(expr{:}); % Get transition matrix from modified txt (controlled network)
    %% Compute Ltilde
    % Needed to convert the system from x(t+1)=Lu(t)x(t) into x(t+1)=Ltx(t)u(t) in order to find u(t) sequence  
    % Mc=Mcontrol(L,n,m); % Get the controllability matrix of the system (calculate only for systems where n<5)
    Lt=L*lwij(2^n,2^m); % Calculate Lt
    %% Get reachability set
    % Depending on initial condition (x0) and final condition (xd) taking from attractors landcape , calculate gene per gene if it possible to reach that finaal state 
    for i1=1:Att_Land_size %Loop 2
        x0=lm(((2^n)-bin2dec(num2str(Att_Land(i1,:)))),2^n);% Initial condition (x0)
        Col_check=(Lt^s)*x0; % According to Theorem 9.3 (Deizhang Cheng et al.2011), the system is reachable from x0 to xd if xd is 
        for i2=1:Att_Land_size %Loop 3
            xd=lm(((2^n)-bin2dec(num2str(Att_Land(i2,:)))),2^n); % Final condition (xd)
            jj=[]; % Matrix jj / It contains all the reachable set 
            uu=1; % Dimensional index of matrix jj
            tt=[]; % Matrix tt / It contains all the reachable xd for each x0 
            kk=1; % Dimensional index of matrix tt
            for nn=1:length(Col_check.v) %Loop 4
                if (Col_check.v(nn)~=jj) % Conditional 1
                    jj(uu,:)=Col_check.v(nn);
                    uu=uu+1; %Increment index uu
                end 
                if (Col_check.v(nn)==xd.v) %Conditional 2
                    ff=lm(nn,2^(m*s));
                    tt(kk,:)=dec2bin(ff.n-ff.v,m*s);
                    kk=kk+1; %Increment index kk
                end
            end 
     %% Get Boolean control sequences 
            if (isempty(tt)==1) %Conditional 3
            else       
                res(zz,:)=[zz i i1 i2]; % Sort information
                U{zz,1}=tt; %Matrix U/ It contains all the possible control sequence for every possible transition
                zz=zz+1; %Increment index zz
            end
        end
    end
    A=regexp(fileread(fullpathname),'\n','split');% Open archive for reading
end
[res2]=SortInfo(res);
[Trajectories_a,Trajectories_b]=ShowTrajectories(res2,Att_Land_size);
assignin('base', 'Trajectories_a', Trajectories_a)
assignin('base', 'Trajectories_b', Trajectories_b)
assignin('base', 'U', U)
set(handles.table_controllability,'Data',Trajectories_a);
set(handles.txt_totalnumbervalue,'String',num2str(length(res2)))
fprintf('Total number of trajectories:\n');
disp(length(res2))
figure
bar(accumarray(res2(:,2), 1))
grid on
xlabel('number of gene');
ylabel('times')
guidata(hObject,handles);


%noise_sig = zeros(size(t));

function varargout = prac5_gui_1(varargin)
% PRAC5_GUI_1 MATLAB code for prac5_gui_1.fig
%      PRAC5_GUI_1, by itself, creates a new PRAC5_GUI_1 or raises the existing
%      singleton*.
%
%      H = PRAC5_GUI_1 returns the handle to a new PRAC5_GUI_1 or the handle to
%      the existing singleton*.
%
%      PRAC5_GUI_1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRAC5_GUI_1.M with the given input arguments.
%
%      PRAC5_GUI_1('Property','Value',...) creates a new PRAC5_GUI_1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before prac5_gui_1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to prac5_gui_1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help prac5_gui_1

% Last Modified by GUIDE v2.5 28-Feb-2022 16:09:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prac5_gui_1_OpeningFcn, ...
                   'gui_OutputFcn',  @prac5_gui_1_OutputFcn, ...
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


% --- Executes just before prac5_gui_1 is made visible.
function prac5_gui_1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prac5_gui_1 (see VARARGIN)

% Choose default command line output for prac5_gui_1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes prac5_gui_1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prac5_gui_1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in dial1.
function dial1_Callback(hObject, eventdata, handles)
% hObject    handle to dial1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "1");
    set(handles.encode, 'string', y);


% --- Executes on button press in dial2.
function dial2_Callback(hObject, eventdata, handles)
% hObject    handle to dial2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "2");
    set(handles.encode, 'string', y);



% --- Executes on button press in dial3.
function dial3_Callback(hObject, eventdata, handles)
% hObject    handle to dial3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "3");
    set(handles.encode, 'string', y);


% --- Executes on button press in dialA.
function dialA_Callback(hObject, eventdata, handles)
% hObject    handle to dialA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "A");
    set(handles.encode, 'string', y);


% --- Executes on button press in dial4.
function dial4_Callback(hObject, eventdata, handles)
% hObject    handle to dial4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "4");
    set(handles.encode, 'string', y);


% --- Executes on button press in dial5.
function dial5_Callback(hObject, eventdata, handles)
% hObject    handle to dial5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "5");
    set(handles.encode, 'string', y);


% --- Executes on button press in dial6.
function dial6_Callback(hObject, eventdata, handles)
% hObject    handle to dial6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "6");
    set(handles.encode, 'string', y);


% --- Executes on button press in dialB.
function dialB_Callback(hObject, eventdata, handles)
% hObject    handle to dialB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "B");
    set(handles.encode, 'string', y);


% --- Executes on button press in dial7.
function dial7_Callback(hObject, eventdata, handles)
% hObject    handle to dial7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "7");
    set(handles.encode, 'string', y);


% --- Executes on button press in dial8.
function dial8_Callback(hObject, eventdata, handles)
% hObject    handle to dial8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "8");
    set(handles.encode, 'string', y);


% --- Executes on button press in dial9.
function dial9_Callback(hObject, eventdata, handles)
% hObject    handle to dial9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "9");
    set(handles.encode, 'string', y);


% --- Executes on button press in dialC.
function dialC_Callback(hObject, eventdata, handles)
% hObject    handle to dialC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "C");
    set(handles.encode, 'string', y);


% --- Executes on button press in dialstar.
function dialstar_Callback(hObject, eventdata, handles)
% hObject    handle to dialstar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "*");
    set(handles.encode, 'string', y);


% --- Executes on button press in dial0.
function dial0_Callback(hObject, eventdata, handles)
% hObject    handle to dial0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "0");
    set(handles.encode, 'string', y);


%--- Executes on button press in dialpound.
function dialpound_Callback(hObject, eventdata, handles)
% hObject    handle to dialpound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "#");
    set(handles.encode, 'string', y);


% --- Executes on button press in dialD.
function dialD_Callback(hObject, eventdata, handles)
% hObject    handle to dialD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    y = get(handles.encode, 'string');
    y = strcat(y, "D");
    set(handles.encode, 'string', y);


% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    set(handles.encode, 'string', "");
    set(handles.decode, 'string', "");
    set(handles.noiseslider, 'value', 0);
    set(handles.noiseamp, 'string', 0);
    set(handles.recsnr, 'string', 0);


% --- Executes on slider movement.
function noiseslider_Callback(hObject, eventdata, handles)
% hObject    handle to noiseslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    %seq = get(handles.encode, 'string');
    f = get(handles.noiseslider, 'value');
    set(handles.noiseamp, 'string', num2str(f));
%     sigamp = 12;
%     noiseamp = f*12;
    %noise_sig = rand(strlen(seq)*size(t));

    


% --- Executes during object creation, after setting all properties.
function noiseslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noiseslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in decodebutton.
function decodebutton_Callback(hObject, eventdata, handles)
% hObject    handle to decodebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    seq = get(handles.encode, 'string');
    % complete this function properly
    Fs = 8e3;
    Lwin = 200;  % can be tweaked
    noise_amp_factor = get(handles.noiseslider, 'value');
    
    fbank = generate_bank(Fs, Lwin);
    [sig, t] = encode(seq, Fs);
    [decoded_seq, snr_out] = decode(sig, t, Fs, fbank, noise_amp_factor);
    set(handles.decode, 'string', decoded_seq);
    set(handles.recsnr, 'string', snr_out);




function h = generate_bank(Fs, Lwin)
    
%Fs = 8e3; % sampling frequency
%t = 0:1/Fs:0.1; % signal duration
%beta = 7e-3;   % can change later
%Lwin = 200; % can be tweaked
    win = 0:1:(Lwin-1);
    
    
    h1 = cos(697*2*pi/Fs*win);% this is my bandpass filter for 697 Hz
    h1 = 1/max(h1) * h1;
    h2 = cos(770*2*pi/Fs*win);
    h2 = 1/max(h2) * h2;
    h3 = cos(852*2*pi/Fs*win);
    h3 = 1/max(h3) * h3;
    h4 = cos(941*2*pi/Fs*win);
    h4 = 1/max(h4) * h4;
    h5 = cos(1209*2*pi/Fs*win);
    h5 = 1/max(h5) * h5;
    h6 = cos(1336*2*pi/Fs*win);
    h6 = 1/max(h6) * h6;
    h7 = cos(1477*2*pi/Fs*win);
    h7 = 1/max(h7) * h7;
    h8 = cos(1633*2*pi/Fs*win);
    h8 = 1/max(h8) * h8;

    h = zeros(8, length(h1));
%     bandfilt = zeros(1, length(filter1));
    
    h(1,:) = h1;
    h(2,:) = h2;
    h(3,:) = h3;
    h(4,:) = h4;
    h(5,:) = h5;
    h(6,:) = h6;
    h(7,:) = h7;
    h(8,:) = h8;
    
%     bandfilt(1,:) = filter1;
%     bandfilt(2,:) = filter2;
%     bandfilt(3,:) = filter3;
%     bandfilt(4,:) = filter4;
%     bandfilt(5,:) = filter5;
%     bandfilt(6,:) = filter6;
%     bandfilt(7,:) = filter7;
%     bandfilt(8,:) = filter8;


function [sig, t] = encode(seq, Fs)            % function for encdoing the input string
    
    m1 = [1 2 3 4 5 6 7 8 9 0 65 66 67 68 42 35];
    m2 = [697 697 697 770 770 770 852 852 852 941 697 770 852 941 941 941];
    m3 = [1209 1336 1477 1209 1336 1477 1209 1336 1477 1336 1633 1633 1633 1633 1209 1477];
    
    M1 = containers.Map(m1, m2);
    M2 = containers.Map(m1, m3);
    
    t = 0:1/Fs:1; %signal time for each segment
    t = 0:1:length(t)-1;    % might be modified later
    
    sig_amp = 12;
    sig = zeros(size(t)); 
    code = zeros(1, strlength(seq));
    
    for i = 1:strlength(seq)
    
%         prompt = '//';
%         x = input(prompt, 's');
        x = seq(i);
        if (x == 'A')
            code(i) = 65;
        elseif (x == 'B')
            code(i) = 66;
        elseif (x == 'C')
            code(i) = 67;
        elseif (x == 'D')
            code(i) = 68;
        elseif (x == '*')
            code(i) = 42;
        elseif (x == '#')
            code(i) = 35;
        else
            code(i) = x - 48;
        end
                    % input has been taken and stored in vector 'seq'
    
        % t = 0:1/Fs:0.1;
        f1 = M1(code(i)); f2 = M2(code(i));
        y = sig_amp*(sin(2*pi*f1/Fs*t) + sin(2*pi*f2/Fs*t));
        
        if(i==1)
            sig = y;
        else
            sig = cat(2, sig, y);
    
        end
    
    end


function [decoded_seq, snr_out] = decode(sig, time, Fs, h, noise_amp_factor)

    len = length(time);
    y = zeros(8, len);
    y_noisefree = zeros(8, len);
    y_omega = zeros(8, len);  % for the DFT's
    
    lookup = [1 2 3 65; 4 5 6 66; 7 8 9 67; 42 0 35 68]; % lookup table    
   
    decoded_seq = "";
    sig_amp = 12; %can be changed accordingly, also for GUI
     
    %noise_amp_factor = get(handles.noiseslider, 'value');
    noise_sig = sig_amp*noise_amp_factor*(randn(size(sig)));
    %snr = 20*log10(rms(sig)/rms(noise_sig));
    temp = sig;
    sig = sig + noise_sig;


    N = length(sig)/len;
    
    for i = 1:N    % computing for each segment of length 'len'
    
        segment = sig((i-1)*length(time)+1 : i*length(time));
        segment_noisefree = temp((i-1)*length(time)+1 : i*length(time));

        y = zeros(8, len);
        y_noisefree = zeros(8, len);

        y_omega = zeros(8, len);
        
            for j = 1:8
                y(j, :) = filtfilt(h(j, :), 1, segment);
                y_noisefree(j, :) = filtfilt(h(j, :), 1, segment_noisefree);
                
            end
            for j = 1:8
                y_omega(j, :) = fftshift(fft(y(j, :), len));
    %             figure
    %             plot(win2, abs(y_omega(j, :))/len);
            end %852 and 1477 should be prominent
    
            y_energy = zeros(2, 4);
    
            rowmax = 0; colmax = 0;
            rowidx = 0; colidx = 0;
            for kappa = 1:4    
                y_energy(1, kappa) = rms(y(kappa,:));
                if(y_energy(1, kappa) > rowmax)
                    rowmax = y_energy(1, kappa);
                    rowidx = kappa;
                end
            end
            
            
            for kappa = 1:4    
                y_energy(2, kappa) = rms(y(4+kappa,:));
                if(y_energy(2, kappa) > colmax)
                    colmax = y_energy(2, kappa);
                    colidx = kappa;
                end
            end
            
            num = lookup(rowidx, colidx);
            if(num == 65)
                decoded_seq = append(decoded_seq, "A");
            elseif(num == 66)
                decoded_seq = append(decoded_seq, "B");
            elseif(num == 67)
                decoded_seq = append(decoded_seq, "C");
            elseif(num == 68)
                decoded_seq = append(decoded_seq, "D");
            elseif(num == 42)
                decoded_seq = append(decoded_seq, "*");
            elseif(num == 35)
                decoded_seq = append(decoded_seq, "#");  
            else
                decoded_seq = append(decoded_seq, num2str(num));
    
            end
    
    end

    snr_out = 20*log10(rms(y)/rms(y - y_noisefree));
   
    
    
    

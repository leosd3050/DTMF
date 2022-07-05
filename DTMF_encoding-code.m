% % EXPERIMENT - 3
% ----------------------------------
clc
close all
prompt = "Enter code length to be transmitted:";
N = input(prompt); % length of input
% Fs = 8e3; % sampling frequency
Ts = 1/Fs; % sampling period
% i = 0;
seq = zeros(1, N);

m1 = [1 2 3 4 5 6 7 8 9 0 65 66 67 68 42 35];
m2 = [697 697 697 770 770 770 852 852 852 941 697 770 852 941 941 941];
m3 = [1209 1336 1477 1209 1336 1477 1209 1336 1477 1336 1633 1633 1633 1633 1209 1477];

M1 = containers.Map(m1, m2);
M2 = containers.Map(m1, m3);

t = 0:1/Fs:1; %signal time for each segment
t = 0:1:length(t)-1;    % might be modified later

% beta = 1;   % can change later
% Lwin = 250; %can be tweaked
% win = 0:1:(Lwin-1);
% 
% h1 = beta*cos(2*pi*697*win); % this is my bandpass filter for 697 Hz
% h2 = beta*cos(2*pi*770*win);
% h3 = beta*cos(2*pi*852*win);
% h4 = beta*cos(2*pi*941*win);
% h5 = beta*cos(2*pi*1209*win);
% h6 = beta*cos(2*pi*1336*win);
% h7 = beta*cos(2*pi*1477*win);
% h8 = beta*cos(2*pi*1633*win);

sig_amp = 12;
sig = zeros(size(t)); 
duty = 10;   %can be tweaked, duty cycle for each key press

y_freq = zeros(N, length(t));
for i = 1:N

    prompt = '//';
    x = input(prompt, 's');
    if (x == 'A')
        seq(i) = 65;
    elseif (x == 'B')
        seq(i) = 66;
    elseif (x == 'C')
        seq(i) = 67;
    elseif (x == 'D')
        seq(i) = 68;
    elseif (x == '*')
        seq(i) = 42;
    elseif (x == '#')
        seq(i) = 35;
    else
        seq(i) = x - 48;
    end
                % input has been taken and stored in vector 'seq'

    % t = 0:1/Fs:0.1;
    f1 = M1(seq(i)); f2 = M2(seq(i));
    y = sig_amp*(sin(2*pi*f1/Fs*t) + sin(2*pi*f2/Fs*t));


    w_ = length(t)*duty/100;

    x_ = rectpuls(t-length(t)/2,w_);
    y = y.*x_;
    y_freq(i, :) = fftshift(fft(y, length(t)));
    if(i==1)
        sig = y;
    else
        sig = cat(2, sig, y);

    end

end

% figure
% plot(t, sig(1, :));
% figure
% plot(t, sig(2, :));
% figure
% plot(t, sig(3, :));
% figure
% plot(t, sig(4, :));
tsig = (0:1:length(sig)-1)/Fs;
% figure
% plot(tsig, sig);
% xlabel('Time (in seconds)');
% ylabel('signal value');

win2 = -Fs/2:Fs/(length(sig)-1):Fs/2;
% for i=1:N
%     subplot(2, 2, i)
%     plot(win2, 2*abs(y_freq(i, :))/length(t));
%     if(i==1)
%         title('press - 1');
%     end
%     if(i==2)
%         title('press - 4');
%     end
%     if(i==3)
%         title('press - 9');
%     end
%     if(i==4)
%         title('press - #');
%     end
%     xlabel('Frequency (in Hz)');
%     ylabel('Scaled freq. magnitude');
% end
figure
plot(win2, 2*abs(fftshift(fft(sig)))/length(t));
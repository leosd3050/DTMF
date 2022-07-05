% % DECODING THE MESSAGE SIGNAL
%clc
close all
len = length(t);
N = length(sig)/length(t);
y = zeros(8, len);
y_omega = zeros(8, len);  % for the DFT's

%thres = 10;    % Bandwidth/2 for bandpass filter

lookup = [1 2 3 65; 4 5 6 66; 7 8 9 67; 42 0 35 68]; % lookup table

%decoded_seq = zeros(1, N);
decoded_seq = "";
alpha = 7.5; % can be tweaked
noise_amp = alpha*sig_amp;   
sig = sig + noise_amp*rand(size(sig));
win3 = -Fs/2:Fs/(length(t)-1):Fs/2;

for i = 1:N    % computing for each segment of length 'len'

    segment = sig((i-1)*length(t)+1 : i*length(t));
    y = zeros(8, len);
    y_omega = zeros(8, len);
    
        for j = 1:8
            y(j, :) = filtfilt(h(j, :), 1, segment);
            
        end
        
        for j = 1:8
            y_omega(j, :) = fftshift(fft(y(j, :), len));
            
%             figure
%             plot(win3, abs(y_omega(j, :))/len);
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
        
%         if(i==1)
        figure
        for j = 1:8
            subplot(2, 4, j)
            plot(win3, abs(y_omega(j, :))/len);
            if(j==1)
                title('fc = 697 Hz');
            end
            if(j==2)
                title('fc = 770 Hz');
            end
            if(j==3)
                title('fc = 852 Hz');
            end
            if(j==4)
                title('fc = 941 Hz');
            end
            if(j==5)
                title('fc = 1209 Hz');
            end
            if(j==6)
                title('fc = 1336 Hz');
            end
            if(j==7)
                title('fc = 1477 Hz');
            end
            if(j==8)
                title('fc = 1633 Hz');
            end
        end
        
        sgtitle(['Decoding pressed key no.: ', num2str(i)]);
        
%         end

end

disp(decoded_seq);

% for i = 1:4
%     figure
%         for j = 1:8
%             subplot(2, 4, j)
%             plot(win3, abs(y_omega(j, :, i))/len);
%             if(j==1)
%                 title('fc = 697 Hz');
%             end
%             if(j==2)
%                 title('fc = 770 Hz');
%             end
%             if(j==3)
%                 title('fc = 852 Hz');
%             end
%             if(j==4)
%                 title('fc = 941 Hz');
%             end
%             if(j==5)
%                 title('fc = 1209 Hz');
%             end
%             if(j==6)
%                 title('fc = 1336 Hz');
%             end
%             if(j==7)
%                 title('fc = 1477 Hz');
%             end
%             if(j==8)
%                 title('fc = 1633 Hz');
%             end
%         end
% 
% end

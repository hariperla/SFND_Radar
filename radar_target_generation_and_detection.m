clear;
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8

%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
tgt_init_range = 50;
tgt_init_velocity = -10;

%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

fc= 77e9; % Operating carrier frequency of Radar 
d_res = 1; % range resolution
speed_of_light = 3e8;
r_max = 200; % max range

% Calculate bandwidth
b_sweep = speed_of_light / (2 * d_res);

% Calculate Tchirp
T_chirp = 5.5 * 2 * r_max / speed_of_light;

% Calculate slope
slope = b_sweep / T_chirp;
                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*T_chirp,Nr*Nd); %total time for samples
L = length(t);

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,L); %transmitted signal
Rx=zeros(1,L); %received signal
Mix = zeros(1,L); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,L);
td=zeros(1,L);


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:L         
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    tgt_range = tgt_init_range + (t(i) * tgt_init_velocity);
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    td = 2 * tgt_range / speed_of_light;
    t_rx = t(i) - td;
    Tx(i) = cos(2 * pi * (fc * t(i) + (0.5 * slope * t(i)^2)));
    Rx (i)  = cos(2 * pi * (fc * t_rx + (0.5 * slope * t_rx^2)));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
    
end

%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix_reshaped = reshape(Mix,[Nr,Nd]);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize. 
fft_Mix_range = fft(Mix_reshaped);

 % *%TODO* :
% Take the absolute value of FFT output
fft_Mix_range = abs(fft_Mix_range/max(max(fft_Mix_range)));

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
fft_one_side = fft_Mix_range(1:(L/2)+1);

%plotting the range
figure ('Name','Range from First FFT')
% subplot(2,1,1)

 % *%TODO* :
 % plot FFT output 
freq = L*(0:L/2)/L;
plot(freq,fft_one_side);
axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix_reshaped=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix_reshaped,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
figure ('Name','2D FFT plot')
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
tc_range = 8;
tc_doppler = 4;

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
gc_range = 4;
gc_doppler = 2;

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 1.6;

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
% noise_level = zeros(1,1);


% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
% CFAR
RDM_thresholded = RDM/max(max(RDM));

for i = tc_range + gc_range + 1 : (Nr/2) - (gc_range + tc_range)
    for j = tc_doppler + gc_doppler + 1 : Nd - (gc_doppler + tc_doppler)
        noise_level = zeros(1,1);
        % Calculate noise sum across the CUT here
        for k = i - (tc_range + gc_range) : i + (tc_range + gc_range)
            for m = j - (tc_doppler + gc_doppler) : j + (tc_doppler + gc_doppler)
                if(abs(i - k) > gc_range || abs(j - m) > gc_doppler)
                    noise_level = noise_level + db2pow(RDM_thresholded(k,m));
                end
            end
        end
        % Determine the training cell size
        training_cell_size = 2 * (tc_doppler + gc_doppler + 1)* 2 * (tc_range + gc_range + 1)-(gc_range * gc_doppler)-1;
        % Add offset to determine the threshold
        threshold = pow2db(noise_level/training_cell_size);
        % Add offset
        threshold = threshold + offset;

        if (RDM_thresholded(i,j) > threshold)
            RDM_thresholded(i,j) = 1;
        else
            RDM_thresholded(i,j) = 0;
        end
    end
end

% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
RDM_thresholded(union(1:(tc_range + gc_range),end - (tc_range + gc_range - 1):end),:) = 0; % setting 0 across rows
RDM_thresholded(:,union(1:(tc_doppler + gc_doppler),end - (tc_doppler + gc_doppler - 1):end)) = 0; % setting 0 across columns
        

% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure ('Name','CFAR output')
surf(doppler_axis,range_axis,RDM_thresholded);
colorbar;

%% Loading audio File and Filter Objects

clc;close all;clear

[y,Fs] = audioread('rec.wav');
load('band_pass_filters.mat');
% load('Filter_1.mat');
% load('Filter_2.mat');
% load('Filter_3.mat');
% load('Filter_4.mat');
% load('Filter_5.mat');
% load('Filter_6.mat');
% load('Filter_7.mat');
% load('Filter_8.mat');
load('Filter_low_400.mat');

Filter_centre=[131 189 272 391 653 810 1165 1676];

%sound(y,Fs);

%% Change sampling rate

if Fs ~= 48000
    [P,Q] = rat(48e3/Fs);
    y = resample(y,P,Q);
    Fs=48000;
end

%% Time axis
length = max(size(y));
maxY = max(abs(y));  % abs(y) is the amplitudes in an all-positive sense

% t = 0:(1/Fs):length/Fs-(1/Fs);

%% considering mono channel and amplify

y = y(:,1)*40;

%% Zero padding 

samples_in_chunk = Fs*0.004;

remainder = mod(length,samples_in_chunk);
padding_length = samples_in_chunk-remainder;

if remainder ~= 0
    padding_zeros = zeros(padding_length,1);
    y(end+1:end+padding_length,1) = padding_zeros;
end

%sound(y,Fs);
%% Original audio Signal plotting

length = max(size(y));
t = 0:(1/Fs):length/Fs-(1/Fs);

limits=length/Fs;

% figure(1)
% plot(t,y(:,1))
% ylabel({'Amplitude'});
% xlabel('Time (s)');
% xlim([0 limits])
% %ylim([-maxY maxY])
% title({'Recorded Signal after padding for Filtering'});

%% design pre-emphasis filter y[n]=y[n-1] - alpha*x[n]

apre_emp=1;
bpre_emp=[1 -0.95];
% fvtool(bpre,apre); % Filter response
ypre_emp=filter(bpre_emp,apre_emp,y);
% figure(2)
% plot(ypre_emp);
% ylabel({'Amplitude'});
% xlabel('Samples');
% xlim([0 length])
% title('Pre-emphasized signal');
%sound(ypre_emp*10,Fs);


%% signal amplification

pre_amp=ypre_emp;

%sound(pre_amp,Fs);

%% Filtering the audio Signal by using Filters 1-8 ( Channel 1 )

Filtered_Signal_From_Filter_1 = filter(Filter_1,pre_amp);
Filtered_Signal_From_Filter_2 = filter(Filter_2,pre_amp);
Filtered_Signal_From_Filter_3 = filter(Filter_3,pre_amp);
Filtered_Signal_From_Filter_4 = filter(Filter_4,pre_amp);
Filtered_Signal_From_Filter_5 = filter(Filter_5,pre_amp);
Filtered_Signal_From_Filter_6 = filter(Filter_6,pre_amp);
Filtered_Signal_From_Filter_7 = filter(Filter_7,pre_amp);
Filtered_Signal_From_Filter_8 = filter(Filter_8,pre_amp);

% plot(t,Filtered_Signal_From_Filter_8)
% ylabel({'Amplitude'});
% xlabel('Time (s)');
% title({'1676 Bandwidth Filter'});

%sound(Filtered_Signal_From_Filter_1 + Filtered_Signal_From_Filter_8,Fs);
%sound(Filtered_Signal_From_Filter_8*10,Fs);
 
%% All filtered signals

Channel(1,:) = Filtered_Signal_From_Filter_1; 
Channel(2,:) = Filtered_Signal_From_Filter_2; 
Channel(3,:) = Filtered_Signal_From_Filter_3; 
Channel(4,:) = Filtered_Signal_From_Filter_4; 
Channel(5,:) = Filtered_Signal_From_Filter_5; 
Channel(6,:) = Filtered_Signal_From_Filter_6; 
Channel(7,:) = Filtered_Signal_From_Filter_7; 
Channel(8,:) = Filtered_Signal_From_Filter_8; 

%sound(Channel(8,:),Fs);

%% Full wave rectification and filter (Envelops) for channels
full_rect_channel=zeros(336768,8);
Filtered_Signal_From_Filter_low_400=zeros(length,8);

for i=1:8
    full_rect_channel(:,i) = abs(Channel(i,:)); 
    Filtered_Signal_From_Filter_low_400(:,i) = filter(Filter_low_400,full_rect_channel(:,i));
end

% figure(1)
% plot(t,full_rect_channel(:,1))
% plot(t,Filtered_Signal_From_Filter_low_400(:,1))

%sound(full_rect_channel(:,1)*10,Fs);
%sound(Filtered_Signal_From_Filter_low_400(:,1)*10,Fs);

%% Quantizing the Filtered signal of the Channel 1 - x_bit
number_of_levels = 2^6;

index_channel=zeros(length,8);
quants_channel=zeros(length,8);

for i = 1:8
    max_val = max(Filtered_Signal_From_Filter_low_400(:,i));
    min_val = 0;
    
    step_size = abs((max_val-min_val)/number_of_levels);
    
    partition = min_val+step_size:step_size:max_val-step_size;
    codebook = 0:step_size:max_val-step_size;
    
    [index_channel(:,i),quants_channel(:,i)] = quantiz(Filtered_Signal_From_Filter_low_400(:,i),partition,codebook);
end

% plot(t,quants_channel(:,1))
% ylabel({'Amplitude'});
% xlabel('Time (s)');
% title({'Quantized Channel-1 signal - 4bit'});


%sound(transpose(quants_channel(:,1))*10,Fs)

%%  rms values of envelop detected signals of channel 

quantized_time_array = 0:1/Fs:0.004-1/Fs;
output_of_channel = zeros(length,8);

for i = 1:8
    for ii = 1:samples_in_chunk:length 
        chunk_4ms = quants_channel(ii:ii+(samples_in_chunk-1),i);
        chunk_rms = rms(chunk_4ms); % RMS calculation
        output_of_channel(ii:ii+(samples_in_chunk-1),i) = sin(2*pi*Filter_centre(i)*quantized_time_array).*chunk_rms; % sin wave assigning
    end
end

% figure(1)
% plot(t,output_of_channel(:,1))
% ylabel({'Amplitude'});
% xlabel('Time (s)');
% title({'Reconstructed Signal'}); 

%sound(output_of_channel(:,1)*10,Fs)


%% Summing all the channels and playing the output 

summed_output = sum(output_of_channel,2);


figure(10)
plot(t,summed_output)
ylabel({'Amplitude'});
xlabel('Time (s)');
xlim([0 limits])
title({'Reconstructed Full  Signal'}); 

%sound(summed_output,Fs)
 
 
 
%% Visualizations

figure(1)
for i=1:8
    subplot(2,4,i) 
    plot(t,transpose(Channel(i,:)))
    ylabel('Amplitude')
    xlabel('Time (s)')
    xlim([0 length/Fs])
    title('Filtered Signal of channel ' + string(i));
end


figure(2)
for i=1:8
    subplot(2,4,i) 
    plot(t,full_rect_channel(:,i))
    ylabel('Amplitude')
    xlabel('Time (s)')
    xlim([0 length/Fs])
    title('Rectified signal of channel ' + string(i)); 
end


figure(3)
for i=1:8
    subplot(2,4,i) 
    plot(t,Filtered_Signal_From_Filter_low_400(:,i))
    ylabel('Amplitude')
    xlabel('Time (s)')
    xlim([0 length/Fs])
    title('Envelope of channel ' + string(i)); 
end

figure(4)
for i=1:8
    subplot(2,4,i) 
    plot(t,quants_channel(:,i))
    ylabel('Amplitude')
    xlabel('Time (s)')
    xlim([0 length/Fs])
    title('Quantized Signal of channel ' + string(i));
end


figure(5)
for i=1:8
    subplot(2,4,i) 
    plot(t,output_of_channel(:,i))
    ylabel('Amplitude')
    xlabel('Time (s)')
    xlim([0 length/Fs])
    title('Reconstructed Full of channel ' + string(i)); 
end



 
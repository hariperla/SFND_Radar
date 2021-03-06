# SFND_Radar
Solution for Radar target generation and detection.

## Initial target range and velocity
```
init_range = 50
init_vel = -10
```

## FMCW Waveform generation
```
c = speed of light
B_sweep = c/(2 * dres)
Tchirp = 5.5 * 2 * max_range / c
slope = B_sweep / Tchirp
```

## Signal generation and moving target simulation
```
td = 2 * tgt_range / speed_of_light;
t_rx = t(i) - td;
Tx(i) = cos(2 * pi * (fc * t(i) + (0.5 * slope * t(i)^2)));
Rx (i)  = cos(2 * pi * (fc * t_rx + (0.5 * slope * t_rx^2)));
```

## Range measurement and Range doppler response
```
Mix_reshaped = reshape(Mix,[Nr,Nd]);
fft_Mix_range = fft(Mix_reshaped);
fft_Mix_range = abs(fft_Mix_range/max(max(fft_Mix_range)));
fft_one_side = fft_Mix_range(1:(L/2)+1);

Mix_reshaped=reshape(Mix,[Nr,Nd]);
sig_fft2 = fft2(Mix_reshaped,Nr,Nd);

sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;
```
## Range measurement output
![image](https://user-images.githubusercontent.com/40438231/159948158-9fe7a9a4-6d94-4c77-9c8a-ec9e94dca6b8.png)


## Range Doppler response output
![image](https://user-images.githubusercontent.com/40438231/159948012-c5720bf7-3277-485e-bc13-f8e7423c34d3.png)


## CFAR implementation
* Slide cell under test across the whole matrix
```
for i = tc_range + gc_range + 1 : (Nr/2) - (gc_range + tc_range)
for j = tc_doppler + gc_doppler + 1 : Nd - (gc_doppler + tc_doppler)
```
* Loop through training cells & guard celsl across range and doppler
```
for k = i - (tc_range + gc_range) : i + (tc_range + gc_range)
for m = j - (tc_doppler + gc_doppler) : j + (tc_doppler + gc_doppler)
```
* Calculate noise level
```
noise_level = noise_level + db2pow(RDM_thresholded(k,m))
```
* Calculate threshold after adding offset
```
threshold = pow2db(noise_level/training_cell_size);
threshold = threshold + offset;
```
* If Cell under test is > threshold then it's 1 else it's 0

## CFAR parameters
```
tc_range = 8
tc_doppler = 4
gc_range = 4
gc_doppler = 2
offset = 2
```

## Steps for suppressing non thresholded cells
```
RDM(union(1:(tc_range + gc_range),end - (tc_range + gc_range - 1):end),:) = 0
RDM(:,union(1:(tc_doppler + gc_doppler),end - (tc_doppler + gc_doppler - 1):end)) = 0
```
## CFAR output
![image](https://user-images.githubusercontent.com/40438231/160125361-7b0dd095-03fa-470b-9219-f569321df6aa.png)


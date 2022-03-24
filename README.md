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

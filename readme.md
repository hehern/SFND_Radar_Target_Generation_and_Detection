# Radar Target Generation and Detection

## FMCW Waveform Generation

```
range_max = 200;
range_reso = 1;
velo_max = 100;
speed_light = 3e8;

fc= 77e9;             %carrier freq
B_sweep = speed_light / (2*range_reso);
T_chirp = 5.5*2*range_max/speed_light;
Slope = B_sweep / T_chirp;
```

## Signal generation and Moving Target simulation

```
for i=1:length(t)         
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = target_position + t(i)*target_velo;
    td(i) = r_t(i)*2/speed_light;
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*(fc*t(i)+0.5*Slope*t(i)^2));
    Rx (i) = cos(2*pi*(fc*(t(i)-td(i))+0.5*(t(i)-td(i))^2*Slope));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
    
end
```

## 2D_FFT

```
Mix = reshape(Mix,[Nr, Nd]);

fft = fft(Mix,Nr);

fft = abs(fft);
fft = fft/max(fft);

fft_out = fft(1:Nr/2+1);

figure ('Name','Range from First FFT')
plot(fft)
axis ([0 200 0 1]);
ylabel('Normalized Amplitude')
xlabel('Range');
```

## RANGE DOPPLER RESPONSE

```
Mix=reshape(Mix,[Nr,Nd]);

sig_fft2 = fft2(Mix,Nr,Nd);

sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);
```

## CFAR implementation

```
Tcr = 10;
Tcd = 5;

Gcr = 4;
Gcd = 2;

offset = 1.6;

noise_level = zeros(Nr/2-2*(Tcd+Gcd),Nd-2*(Tcr+Gcr));
number = (2*(Tcd+Gcd)+1)*(2*(Tcr+Gcr)+1)-((2*Gcr+1)*(2*Gcd+1));

CFAR= zeros(size(RDM));
for j=1:Nd-2*(Tcr+Gcr)
    for i=1:Nr/2-2*(Tcd+Gcd)      
        tcp = db2pow(RDM(i:i+2*(Tcd+Gcd),j:j+2*(Gcr+Tcr)));
        tcp(Tcd+1:end-Tcd,Tcr+1:end-Tcr) = 0;
        
        noise_level(i,j) = pow2db(sum(sum(tcp))/number);
        thresh = noise_level(i,j)*offset;
        if RDM(i+(Tcd+Gcd),j+(Tcd+Gcr))>thresh
            CFAR(i+(Tcd+Gcd),j+(Tcd+Gcr)) = 1;
        else
            CFAR(i+(Tcd+Gcd),j+(Tcd+Gcr)) = 0;
        end
           
    end
end
```
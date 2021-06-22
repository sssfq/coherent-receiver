% optical system data extraction
file = "PDM16QAM_28GBd_3600km_-3dBm_17.9dB_OSNR.bin";
f = fopen(file);
line1 = fgets(f);
line2 = fgets(f);
line3 = fgets(f);
fclose(f);
[par,val] = strtok(line3, ':');
lines = eval(val(3:end));
f = fopen(file);
header = [];
CurrMainField = '';  
for i=1:lines
    s = fgets(f);
    if ( s(1)=='[')
    % We found section
        CurrMainField = genvarname(lower(s(2:end-2)));
        header.(CurrMainField) = [];    % Create field in Result
    else
        [par,val] = strtok(s, ':');
        val = val(3:end-1);
        if ~isempty(CurrMainField)
            % But we found section before and have to fill it
            header.(CurrMainField).(lower(genvarname(par))) = val;
        else
            % No sections found before. Orphan value
            header.(lower(genvarname(par))) = val;
        end
    end
end
m = eval(header.fileformat.data_ncolumns);
n = eval(header.fileformat.data_nrows);
data = fread(f,[m n],eval(header.fileformat.data_precision));
fclose(f);

Xi = data(1,:);
Xq = data(2,:);
Yi = data(3,:);
Yq = data(4,:);

InXo = Xi + 1i * Xq;
InYo = Yi + 1i * Yq;
InXo = InXo-mean(InXo);
InYo = InYo-mean(InYo);

SymbolRate = 28e9;
SampleRate = 50e9;
SamplesPerSymbol = SampleRate/SymbolRate;
Tsam = 1/SampleRate; %采样间隔
Tsaml = 1/(SymbolRate*2);%二分之一倍的码元间隔
InXs = interp1(0:Tsam:(length(InXo)-1)*Tsam, InXo, 0:Tsaml:(length(InXo)-1)*Tsam,'spline');
InYs = interp1(0:Tsam:(length(InXo)-1)*Tsam, InYo, 0:Tsaml:(length(InXo)-1)*Tsam,'spline');
c0 = 299792458;
ChannelWavelength = 1550.92 * 1e-9;
L = 3600e3;
S = 0.07 * 1e-12 / 1e-15;
D = -16.75 * 1e-12 / 1e-6;

NumD_Eff = 9999;                                                         % Odd number
DTime    = Tsaml;%二分之一倍的码元间隔
DFreq    = (1 / DTime) / (NumD_Eff - 1);%频率间隔=两倍码元频率等分NumD_Eff-1份
VFreq    = (-(NumD_Eff-1) / 2 : (NumD_Eff-1) / 2) * DFreq;%-的码元频率到+的码元频率 序列
VOmeg    = 2 * pi * VFreq;                                               % Angular frequency 

Hf       = exp(-1i * (VOmeg .^ 2) * ( (ChannelWavelength ^ 2 * D * L) / (4 * pi * c0) ));% + (1i * ( (S * L * ChannelWavelength ^ 4 * VOmeg .^ 3) / (24 * pi ^ 2 * c0 ^ 2) )) );
Ht       = (fftshift(ifft(ifftshift((Hf)))));
Htc      = conj(Ht);

B = SymbolRate * 2; %带宽 

TapNumberD = 2 * floor((abs(D) * (ChannelWavelength ^ 2) * L) / (2 * c0 * Tsaml^2)) + 1; %计算抽头数
if mod(TapNumberD,2) == 0
    TapNumberD = TapNumberD + 1;
end

TapD        = Htc(1,((NumD_Eff+1)/2-(TapNumberD-1)/2):((NumD_Eff+1)/2+(TapNumberD-1)/2));

OutX_Temp   = conv(InXs,TapD);           
conL = length(OutX_Temp);

OutX        = OutX_Temp(1,((TapNumberD-1)/2) + 1 : end - ((TapNumberD-1)/2));
papr(OutX_Temp)
OutY_Temp   = conv(InYs,TapD);           
OutY        = OutY_Temp(1,((TapNumberD-1)/2) + 1 : end - ((TapNumberD-1)/2));


figure(1)
fftx = fft(OutX);
ffty = fft(OutY);
xrfx = abs(xcorr(fftx));
yrfx = abs(xcorr(fftx,ffty));
plot(yrfx(1:(end-1)/2)+xrfx(1:(end-1)/2))

figure(2)
fftxi=fft(InXo);
fftyi=fft(InYo);
xirfx = abs(xcorr(fftxi));
yirfx = abs(xcorr(fftxi,fftyi));
plot(yirfx(1:(end-1)/2)+xirfx(1:(end-1)/2))

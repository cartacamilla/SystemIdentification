function model = fourier_analysis(y,u,Te,average)

N = length(u);

xcorr_uu = intcor(u,u);
[pks,locs]=findpeaks(xcorr_uu, 'MinPeakheight',0.05);
u_period = mean(diff(locs));
nb_period = N / u_period;

%Fourier transform
w_s = 2*pi/Te;
w_nyquist = w_s/2;

NYQUIST_INDEX = round(u_period/2);

if average == true
    avg_y = zeros(u_period,1);
    for i = 1:nb_period
        i_0 = 1+(i-1)*u_period;
        i_1 = i*u_period;
        
        avg_y = avg_y + fft(y(i_0:i_1));
    end
    
    fftY = avg_y ./ nb_period;
    
    %no average for u beacuse it's an open loop system
    fftU = fft(u(1:u_period));
    
    
    g = fftY./fftU;
    Ng = length(g(1:end/2));
    w_n = w_nyquist*(0:(Ng-1))/Ng;

    model = frd(g(1:Ng), w_n);
    
else
    fftY = fft(y);
    
    fftU = fft(u);
    
    g = fftY./fftU;
    Ng = length(g(1:end/2));
    
    w_n = w_nyquist*(0:(Ng-1))/Ng;

    model = frd(g(1:Ng), w_n);    
end

h1 = figure(1);
bode(model, w_n)

title('Bode Diagram');
legend('Identified model - Fourier analysis with averaging');
saveas(h1, '../images/1_Fourier_analysis', 'png');

h2 = figure(2);

Y1 = fft(y(1:u_period));
U1 = fft(u(1:u_period));
Ny = length(Y1);
w_n = w_nyquist*(0:(Ny-1))/Ny;

Y1_model = frd(Y1(1:end/2), w_n(1:2:end-1));
U1_model = frd(U1(1:end/2), w_n(1:2:end-1));

bode(Y1_model,U1_model, w_n);

title('Bode Diagram');
legend('Fft of output y', 'Fft of input u');
saveas(h2, '../images/2_Fourier_analysis', 'png');



end
function model = spectral_analysis(y,u,Te,m,average,window_type)

N = length(u);

if nargin < 4
    average = false;
end

if nargin < 5
    window = ones(N,1);
end

size_split = floor(N/m);

M = 550;

switch window_type
    case 'hann'
        window = hann(M);
    case 'hamming'
        window = hamming(M);
end

%correlation
PHI_yu = zeros(2*size_split,1);
PHI_uu = zeros(2*size_split,1);

window_pad = [zeros(ceil((2*size_split-M-1)/2),1);window;zeros(floor((2*size_split-M)/2),1)];

if average == true
    for i = 1:m
        i_0 = 1+(i-1)*size_split;
        i_1 = (i)*size_split;

        PHI_yu = PHI_yu + fft(window_pad .* intcor(u(i_0:i_1), y(i_0:i_1)));
        PHI_uu = PHI_uu + fft(window_pad .* intcor(u(i_0:i_1), u(i_0:i_1)));

    end
    PHI_yu = PHI_yu./m;
    PHI_uu = PHI_uu./m;

    G = PHI_yu./PHI_uu;

    Ng = length(G(1:end/2));

    w_s = 2*pi/Te;
    w_nyquist = w_s/2;
    w_n = w_nyquist*(0:Ng-1)/Ng;

    model = frd(G(1:Ng), w_n);



end

h3 = figure(3);

bode(model,w_n);

title('Bode diagram')
legend('Identified model - Spectral analysis with averaging and windowing')
saveas(h3, 'final-report/images/1_Spectral_analysis', 'png');


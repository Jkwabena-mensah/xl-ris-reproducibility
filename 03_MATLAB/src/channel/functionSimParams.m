function p = functionSimParams()
%Return a struct with all shared simulation parameters for the
%PW-BCD-MM near-field XL-RIS beam-tracking paper (Table II), extended with
%the blocked-direct-link + LoS-dominant NLOS channel parameters.
%
%This is version 2.1 (Last edited: 2026-06-02)  -- LoS-dominant (kappa=15 dB)
%License: GPLv2.

%%Physical / array parameters
p.c       = 3e8;                 %Speed of light [m/s]
p.fc      = 10e9;                %Carrier frequency [Hz]
p.lambda  = p.c / p.fc;          %Wavelength [m] (= 0.03 m)
p.N1      = 16;                  %RIS rows
p.N2      = 16;                  %RIS columns
p.N       = p.N1 * p.N2;         %Total RIS elements (= 256)
p.D       = 1.20;                %RIS aperture (largest dimension) [m]
p.dR      = 2 * p.D^2 / p.lambda;%Rayleigh distance [m] (= 96 m)
p.dElem   = p.D / (p.N1 - 1);    %Inter-element spacing [m]

%%Network geometry
p.M       = 4;                   %BS antennas (absorbed into effective feed)
p.K       = 8;                   %Number of single-antenna mobile users
p.pBS     = [30; 0; 5];          %BS position [m]
p.region  = [3 44; -5 5; 0 10];  %Service region R [m]

%%Link budget
p.SNRdB   = 0;                   %Operating SNR [dB]
p.sigma2  = 1;                   %Normalised noise variance
p.Ptx     = 10^(p.SNRdB/10) * p.sigma2;  %Transmit power (SNR = Ptx/sigma2)
p.alpha   = 1 / p.N;             %Channel amplitude alpha_k = 1/N

%%Codebook (region-gridded near-field focal points)
p.cbNx    = 16; p.cbNy = 10; p.cbNz = 16;   %16*10*16 = 2560 codewords

%%BCD-MM solver settings
p.epsTol  = 1e-4;                %Outer-loop tolerance
p.muTol   = 1e-2;                %Inner MM tolerance (calibrated 2026-07-23: former
                                 %1e-5 was ~6e-7 rad/element, far tighter than needed
                                 %and forced the inner loop to hit Imax=100 every time,
                                 %so nIter pinned at the ~5000 cap. At 1e-2 (0.036 deg
                                 %precision) the rate is identical but the solver
                                 %converges in ~hundreds of updates, giving a clean,
                                 %cap-independent ~50%% warm-start iteration saving.)
p.Lmax    = 50;                  %Max outer iterations
p.Imax    = 100;                 %Max inner MM iterations

%%Mobility / slot
p.Ts      = 0.5e-3;              %Slot duration [s]
p.maxTurn = 30*pi/180;           %Max heading change per slot [rad]

%%EKF
p.q       = 1e-4;                %Process-noise variance [m^2]
p.r       = 1e-3;                %Measurement-noise variance [m^2]

%%Near-field focal beamwidth (for normalising prediction error)
p.zNom    = 25;                  %Nominal focal depth [m]
p.beamwidth = p.lambda * p.zNom / p.D;   %-3 dB focal width ~ lambda*z/D

%%Extended channel: blocked direct link + LoS-dominant NLOS multipath
p.kappaDB   = 15;                  %Near-field Rician factor [dB] (LoS-dominant)
p.kappa     = 10^(p.kappaDB/10);   %... in linear scale for the channel model
p.Lscat     = 8;                   %NLOS scatterers per user
p.scatSpread= 2.0;                 %Std of scatterer offset around a user [m]
p.rhoDirect = 0;                   %Direct BS->UE amplitude (0 = blocked by obstacle)

end

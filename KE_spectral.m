function [KE_operator, KE_energy] = KE_spectral(SPACE, TIME, Psi_R, Psi_I)
%% Desc.
%**************************************************************************
% Kinetic Operator
%--------------------------------------------------------------------------
% Pointwise multiplication:
% 1)    FFT --> Fourier Space
% 
% 2)    Plane Wave expansion in fourier domain
% 
% 3)    Inverse FFT --> Position Space
%
% for vector potential, just use "(Kvec + A_pot).^2/2" for wave vector
% T_op = diag(Kwave);
%*(hPlanck^2/me)
%__________________________________________________________________________



%% Sanity Check





%% Local Defs.
CONSTANTS;







%% Kinetic Energy from Fourier Domain
Kvec = fftaxisshift(fftaxis(SPACE.Axis));
Kwave = (ifft(2*pi*((Kvec).*fft(Psi_R+1i*Psi_I)))).^2*(hPlanck^2/(2*me));
% (hPlanck^2/(2*me))

% 2*pi*(ifft((Kvec.^2*(-hPlanck^2)/(2*me)).*fft(Psi_R+1i*Psi_I)));


%% Output
% KE_real = real(Kwave);
% KE_imag = imag(Kwave);

KE_operator = Kwave;
KE_energy = trapz(Kwave,SPACE.Axis);



end
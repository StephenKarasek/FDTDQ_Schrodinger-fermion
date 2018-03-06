function [KE_operator, KE_energy] = KE_gradient(SPACE, TIME, Psi_R, Psi_I)
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







%% Kinetic Energy Operator



Kin_ = (-(hPlanck^2)/(2*me))*...
    sum((Psi_R-1i*Psi_I).*...
    del2(Psi_R,SPACE.Axis) + ...
    del2(1i*Psi_I,SPACE.Axis));

Kin = (-(hPlanck^2)/(2*me))*...
    ((Psi_R-1i*Psi_I).*...
    del2(Psi_R,SPACE.Axis) + ...
    del2(1i*Psi_I,SPACE.Axis));


% Kin = gradient(gradient(Psi_R+1i*Psi_I,SPACE.delta),SPACE.delta).*(-hPlanck^2)/(2*me);


% Remove Artefacts & Spurrious Reflections
% Knorm = sum(Kin)/numel(Kin);
% for n=1:SPACE.N; if (Kin(n)>1e1*Knorm); Kin(n) = Knorm; end; end



%% Output
% KE_real = real(Kin);
% KE_imag = imag(Kin);

KE_operator = Kin;

% KE_energy = trapz(SPACE.Axis_, Kin);
% KE_energy = trapz(SPACE.Axis, Kin);
KE_energy = Kin_;

end
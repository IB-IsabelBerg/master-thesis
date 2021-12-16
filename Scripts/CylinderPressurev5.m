function[P_a] = CylinderPressurev5(ka,M,phi_vec)
% P_a = vector of calculated pressure values
% k = wave number: (2*pi*f)/c, 
% f = frequency of plane incoming wave [Hz]
% c = velocity of sound [m/s]
% a = radius of cylinder [m]
% M = number of iterations in sum
% phi_vec = vector of azimuth angles that the pressure is calculated for

% Function that calculates the pressure on the surface of an ideal,
% infinitely long cylinder. Equation (11.29) from Junger and Feit *
% Written April 8th 2021 by Isabel Berg

%% Calculations
% Pressure amplitude can be set equal 1, since the eventual division in the
% directivity factor removes this factor, since it's the same in both the
% numerator and the denominator.
P0 = 1;%sqrt(rho*c*I); %Pressure amplitude of incoming plane wave
epsilon_0 = 1;
epsilon_m = 2;

m_sum = zeros(length(phi_vec),M+1);
for m = 0:M
    H_n_prime = ((m)*besselh((m),1,ka))/(ka)-besselh((m+1),1,ka);
    if m == 0
        b_vec = ((epsilon_0*((1i)^(m+1)))/H_n_prime)*cos(m*phi_vec);
    else
        b_vec = ((epsilon_m*((1i)^(m+1)))/H_n_prime)*cos(m*phi_vec);
    end
    m_sum(:,(m+1)) = b_vec;
end

X = cumsum(m_sum,2); %Summing over M iterations
X_final = X(:,end);

P_a = ((2*P0)/(pi*ka))*X_final;

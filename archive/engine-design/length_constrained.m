% ASSUME CENTRE-PERFORATED FUEL GRAIN %

n = 0.5; % regression rate exponent, empirical, Ch.6, pp. 110
PDR = 2:0.01:3.5; % port dia ratio Df/Di, Ch. 7, pp. 134
D_f = 60e-3; % final fuel grain port dia, m
D_i = D_f./PDR; % initial fuel grain port dia, m
L_f = (50:1:200)*1e-3; % fuel length, m
m_f = zeros(100,100); % total fuel mass, kg
p_f = 900; % fuel density, kg m-3
m_o = 1.5; % total oxidiser mass, kg - GUESS
a_o = 1.55E-04; %regression rate coefficient (only in terms of oxidiser mass flux), empirical, Ch.6, pp. 110

for i = 1:length(PDR)  % Loop over each initial diameter D_i(i)
    for j = 1:length(PDR) % Loop over each fuel length L_f(j)
        m_f(i, j) = 1/4 * (L_f(j) * (pi * p_f * (D_f^2 - D_i(i)^2)));
        
    end
end

%m_f = 1/4 * (L_f .* (pi * p_f * (D_f^2 - D_i.^2)));

t_b = (((D_f^(2*n+1) - D_i.^(2*n+1)) .* pi^n) ./ (a_o * (2*n + 1) * 2^(2*n + 1) * m_o .^ n)) .^ (1 / (1-n));
 

s = surf(L_f, t_b, m_o./m_f, 'FaceAlpha',0.5);

min_OF = 6.7;
max_OF = 9.1;

clim([min_OF max_OF]);
cmap = parula(256); 

cmap(1, :) = [1 0 0]; 
cmap(end, :) = [1 0 0]; 

colormap(cmap);

xlabel('Fuel Grain Length / m'); ylabel('Burn Time / s'); zlabel('O/F Ratio'); colorbar; s.EdgeColor = 'none';
% ASSUME CENTRE-PERFORATED FUEL GRAIN %

n = 0.5; % regression rate exponent, empirical, Ch.6, pp. 110
D_i = 0; % initial fuel grain port dia, m
D_f = 60e-3; % final fuel grain port dia, m
L_f = 0; % fuel length, m
m_f = 0.2; % total fuel mass, kg
p_f = 900; % fuel density, kg m-3
m_o = 0.5; % total oxidiser mass, kg - GUESS
G_of = 0.5; % end of burn mass flux, kg s-1 - GUESS
a_o = 1.55E-04; %regression rate coefficient (only in terms of oxidiser mass flux), empirical, Ch.6, pp. 110

PDR = 2:0.1:3; % port dia ratio Df/Di, Ch. 7, pp. 134
D_i = D_f ./ PDR;
%(((2*n + 1) * 8*a_o * m_o * G_of^(n-1)) / (pi * (1 - (PDR(i)^(-1))^(2*n+1))))^(1/3);

L_f = (4*m_f)./(pi * p_f * (D_f^2 - D_i.^2));



%m_dot_o = (pi * G_of * D_f^2)/4;

%t_b = m_o/m_dot_o;


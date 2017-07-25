function COEFFICIENTS = Initialize_parameters(Number_of_Nucleotides,jitter)
    r_crRNA = 10*(1+jitter*rand(1)*(-1)^(round(1/rand(1))));
    Delta_crRNA = 1*(1+jitter*rand(1)*(-1)^(round(1/rand(1))));
    r_Cas9 = 10*(1+jitter*rand(1)*(-1)^(round(1/rand(1))));
    Delta_Cas9 = .1*(1+jitter*rand(1)*(-1)^(round(1/rand(1))));
    Delta_crRNACas9 = .1*(1+jitter*rand(1)*(-1)^(round(1/rand(1))));
    Lambda =  100*(1+jitter*rand(1)*(-1)^(round(1/rand(1))));
    k_f = 15*(1+jitter*rand(1)*(-1)^(round(1/rand(1))));
    k_I = 15*(1+jitter*rand(1)*(-1)^(round(1/rand(1))));
    D = 250*(1+jitter*rand(1)*(-1)^(round(1/rand(1))));
    V = 1*(1+jitter*rand(1)*(-1)^(round(1/rand(1))));
    kd = 5*(1+jitter*rand(1)*(-1)^(round(1/rand(1))));
    k_c = .0016*(1+jitter*rand(1)*(-1)^(round(1/rand(1))));
    mu = .4*(1+jitter*rand(1)*(-1)^(round(1/rand(1))));
    %initial ODE system coefficients
    COEFFICIENTS1 = [r_crRNA Delta_crRNA k_f r_Cas9 Delta_Cas9 Delta_crRNACas9 k_I D Lambda V kd k_c mu];
    %initial gamma coefficients
    COEFFICIENTS2 = (zeros(1,Number_of_Nucleotides)+.005.*ones(1,Number_of_Nucleotides).*((1+jitter.*rand(1,Number_of_Nucleotides).*(-1).^(round(1./rand(1,Number_of_Nucleotides))))));
    COEFFICIENTS2 = flip(sort(COEFFICIENTS2));
    % binding energy coefficients (both match and mismatch)
    COEFFICIENTS3 = ones(1,10).*((1+jitter.*rand(1,10).*(-1).^(round(1./rand(1,10)))))+1;
    %Initialize Dimention scaling
    COEFFICIENTS4 = 2.*((1+jitter.*rand(1,1).*(-1).^(round(1./rand(1,1)))));
    COEFFICIENTS5 = 2.*((1+jitter.*rand(1,1).*(-1).^(round(1./rand(1,1)))));
    %%%%======================================================================
    %Initialize Parameters
    %%%%======================================================================
    ODE_PARAMETERS = COEFFICIENTS1;
    GAMMA_PARAMETERS = COEFFICIENTS2;
    MATCH_PARAMETERS = COEFFICIENTS3;
    DIMENSION = COEFFICIENTS4;   
    ANCHOR = COEFFICIENTS5;
    COEFFICIENTS = [ODE_PARAMETERS GAMMA_PARAMETERS MATCH_PARAMETERS DIMENSION ANCHOR];
end
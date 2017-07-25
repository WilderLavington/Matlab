function cost = Obj_d(dim)
    global RELATIVE_REPRESSION 
    global SG_RNA 
    global HOST_SITE
    global TEST_NUMBER
    global global_COEFFICIENTS
    %%%%%========================================
    %Initialize parameters
    %%%%%========================================
    global Number_of_Nucleotides
    ODE_Parameters = global_COEFFICIENTS(1:13);
    drop_out = global_COEFFICIENTS(14:14+Number_of_Nucleotides-1)';
    match_mismatch_params = global_COEFFICIENTS(14+Number_of_Nucleotides:14+Number_of_Nucleotides+10-1);
    anchor = global_COEFFICIENTS(14+Number_of_Nucleotides+10);
    dimention = dim;
    tests = unique(TEST_NUMBER);
    cost = 0;
    %%%%%========================================
    %Iterate through all experiments
    %%%%%========================================
    for kk = 1:length(tests)
        test_batch = find(TEST_NUMBER == tests(kk));
        batch_size = length(test_batch);
        P = zeros(batch_size,1);
        for ii = 1:batch_size
            %%%%%========================================
            %calculate binding probability
            %%%%%========================================
            N_genome = [HOST_SITE(test_batch(ii),:); SG_RNA(test_batch(ii),:)]';
            [s,r,l,d] = model_probabilities(N_genome,drop_out,match_mismatch_params,dimention,anchor);
            P(ii) = system_solve_mc(l,r,s,d);
        end
        %%%%%========================================
        %Use stiff, low order numerical solver ODE 23s
        %%%%%========================================
        scaling = 1;
        IC1_to_4 = [scaling,scaling,0,0];
        IC5 =  ones(size(test_batch'));
        Initial_conditions = [IC1_to_4 IC5];
        %interval  = [0,1200];
        tspan = linspace(0,100,500);
        %%%%%========================================
        %create globals for function call in ODE solver
        %%%%%========================================
        global ODE_Parameters_current
        ODE_Parameters_current = ODE_Parameters;
        global Probability_of_binding
        Probability_of_binding = P;
        global N_total_current
        N_total_i = IC5;
        N_total_current = N_total_i;
        f = @Stiff_Solve;
        [t,N5] = ode23s(f,tspan,Initial_conditions');
        new_cleavage = ones(size(N5(:,5:end)));
        for ii = 1:max(size(N5))
            for jj = 5:min(size(N5))
            new_cleavage(ii,jj-4) = (N_total_current(jj-4)-N5(ii,jj))/N_total_current(jj-4);
            end
        end
        complete_repression = new_cleavage(end,1:end)';
        %%%%%========================================
        %calculate cost 
        %%%%%========================================
        cost = cost + sum((RELATIVE_REPRESSION([test_batch])-complete_repression).^2);
    end
end

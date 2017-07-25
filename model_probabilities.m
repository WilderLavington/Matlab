function [s,r,l,d] = model_probabilities(N_genome_i,drop_out_probability,match_mismatch_params,dimention,anchor)
    %%%========================================================================
    %%%========================================================================
    %find what type of binding event each match/mismatch is in sequence then
    %pull the correct coefficients from the match_mismatch_params 
    %%%========================================================================
    %%%========================================================================
    %required keys
    global Number_of_Nucleotides
    sequence_key = [0 0; 0 1; 0 2; 0 3; 1 0; 1 1; 1 2; 1 3; 2 0; 2 1; 2 2; 2 3; 3 0; 3 1; 3 2; 3 3];
    reduced_key = [1 2 3 4 2 5 6 7 3 6 8 9 4 7 9 10; 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
    param_key = match_mismatch_params;
    sphere_radius = zeros(Number_of_Nucleotides,1);
    %assign sphere radius
    mismatches = mismatch_locations(N_genome_i'); % 1 = match, 0 = mismatch
    for ii = 1:Number_of_Nucleotides
        for jj = 1:16
            if  N_genome_i(ii,1) == sequence_key(jj,1) &&  N_genome_i(ii,2) == sequence_key(jj,2)
                sphere_radius(ii) = param_key(reduced_key(1,jj));
            end
        end
    end
    % padding ~ 
    sphere_radius = [1; sphere_radius; anchor];
    mismatches = [1 mismatches' 1];
    drop_out_probability = [1 drop_out_probability' 0];
    %%%========================================================================
    %%%========================================================================
    %calculate cumulative gamma probability given from mismatches
    %%%========================================================================
    %%%========================================================================
    gamma_mismatches = -1.*(mismatches-1);
    %create cumulative gamma matrix
    cum_g_mat = zeros(length(mismatches),length(mismatches));
    for ii = 1:length(mismatches)
        %set relevent values 
        for kk = 1:ii
            cum_g_mat(ii,kk) = 1;
        end
        %mask with positions that have mismatches
        cum_g_mat(ii,:) = gamma_mismatches.*cum_g_mat(ii,:);
    end
    %use cumulative matrix to calculate cumulative instability
    state_gamma = zeros(length(nonzeros(mismatches))+2,1);
    for ii = 1:length(gamma_mismatches)-1
        state_gamma(ii+1) = sum(drop_out_probability.*cum_g_mat(ii,:));
    end
    %%%========================================================================
    %%%========================================================================
    %calculate transition probabilities for model
    %%%========================================================================
    %%%========================================================================
    %inialize vectors
    binding_l = zeros(length(nonzeros(mismatches))+2,1);
    unbinding_r = zeros(length(nonzeros(mismatches))+2,1);
    nomovement_s = zeros(length(nonzeros(mismatches))+2,1);
    %solve for probabilities
    %%%========================================================================
    %%%========================================================================
    %calculate transition probability
    %%%========================================================================
    %%%========================================================================\
    % initializations
    current_state = 0;
    ii = 1;
    current_start = ii;
    % interior
    while ii < length(mismatches)
        %interior mismatches
        if mismatches(ii+1) == 0
            ii = ii + 1;
        %next state
        elseif mismatches(ii+1) == 1 
            ii = ii + 1;
            current_end = ii;
            current_state = current_state + 1;
            %assign probabilities
            binding_l(current_state) = (1-state_gamma(ii))*(sphere_radius(current_end)^dimention)...
                                    /((2*sum(sphere_radius(current_start+1:current_end)))^dimention);                              
            %make sure there is a mismatch between the two sites
            if current_end - current_start >= 1
                unbinding_r(current_state)  = (1-state_gamma(ii))*(sum(sphere_radius(current_start+1:current_end-1).^dimention))...
                                    /((2*sum(sphere_radius(current_start+1:current_end)))^dimention);
            else
                unbinding_r(current_state) = 0;
            end
            nomovement_s(current_state)  = (1-state_gamma(ii))*(1-binding_l(current_state)-unbinding_r(current_state));
            current_start = ii;
        end
    end
    % end nucleotide
    r = [0; unbinding_r(1:current_state); 0];
    l = [0; binding_l(1:current_state); 0];
    s = [1; nomovement_s(1:current_state); 1];
    d = state_gamma;
end
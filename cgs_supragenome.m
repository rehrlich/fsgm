function cgs_supragenome(name, counts)
%Calculate various parameters using the finite supragenome model of Hogg et al.
%The function ''cgs_supragenome'' uses the observed frequencies of genes in a
%set of bacterial genomes and the finite supragenome model of Hogg et al.
%(Genome Biol. 8, R103 (2007)) to calculate the maximum-likelihood (ML) values
%for N, the size (in genes) of a bacterial supragenome (pan-genome); pi, the
%(K-element) vector of mixture coefficients or probabilities that a given gene
%in N will fall in one of the K population gene-frequency-defined classes of genes;
%and mu, the (K-element) vector of population gene-frequencies for each gene class.
%Other statistics are calculated based on the ML values of N, pi, and mu.  One mode
%of running this program is to just calculate these statistics based on previously
%determined ML values of N, pi, and mu.

programStart = datestr(now);
warning('off', 'MATLAB:nearlySingularMatrix')

%{
**************************************************************************
SECTION 1.  GET INPUT DATA FROM USER
**************************************************************************
%}
get_prompts();

fileName_diary = ['CommandWindow_', name, '.txt'];
diary(fileName_diary);

gene_freq_histogram_predicted_strains = size(counts, 2);
new_core_genes_predicted_strains = size(counts, 2) * 10;

fileName_N_vs_likelihood = ['N_vs_likelihood_', name, '.txt'];

gene_freq_histogram_observed = counts;
gclass_probs_start = [0.1 0.3 0.5 0.7 0.9 1.0];
gclass_probs_separation_min = 0.05;
gclass_probs_min = 0.01;
gclass_mixcoeff_min =0.01;
N_start_offset = 10;
N_range = 4000;

%Echo input data (notice the absence of semi-colon line endings):
fileName_N_vs_likelihood
gene_freq_histogram_observed
gclass_probs_start
gclass_probs_separation_min
gclass_probs_min
gclass_mixcoeff_min
N_start_offset
N_range
gene_freq_histogram_predicted_strains
new_core_genes_predicted_strains

%{
**************************************************************************
SECTION 2.  RUN OPTIMIZATION
**************************************************************************
%}

% initialize variables
% number of genomes
S = length(gene_freq_histogram_observed) - 1;
% number of gene frequency classes
K = length(gclass_probs_start);

% calculate a constant factored out of ML optimization funtion:
fixFact = 0;
for i = 1 : S
    for n = 2:gene_freq_histogram_observed(i+1)
        fixFact = fixFact - log(n);
    end
end

% set minimum and maximum allowable value for N (supragenome size)
Nmin = sum( gene_freq_histogram_observed(2:end) ) + N_start_offset;
Nmax = Nmin + N_range;

% starting values for mixing coefficients and gene class probabilities
start_params = [ ones(K,1) ./ K; gclass_probs_start.' ];

% set options for minimization algorithm
options = optimset('LargeScale', 'off', 'FunValCheck', 'on', 'Display', 'off');

% define parameter constraints

% sum of mixture coefficients = 1
% set 'core' gclass_prob = 1
Aeq = [ ones(1,K),   zeros(1,K);      ...
    zeros(1,K),  zeros(1,K-1), 1   ];
beq = [ 1; 1 ];

% make sure Freq[class(i)] < Freq[class(i+1)]
A = [  zeros(K-1,K), diag( ones(K-1,1) ) + diag( -ones(K-2,1), 1), [zeros(K-2,1); -1]  ];
b = -gclass_probs_separation_min * ones(K-1,1);

% bounds mixing coefficients and class frequencies
lb = [ gclass_mixcoeff_min*ones(K,1); gclass_probs_min*ones(K-1,1); 1 ];
ub = [ ones(K,1);                      ones(K,1)                     ];

%Echo data (notice the absence of semi-colon line endings):
Aeq
beq
A
b
lb
ub

%% Optimization

X = [];
ML_logProb = -Inf;
ML_params = start_params;

% begin optimization at N=Nmin
for N = Nmin : Nmax
    % unseen
    gene_freq_histogram_observed(1) = N - ...
        sum( gene_freq_histogram_observed( 2:length(gene_freq_histogram_observed) ) );
    fact = fixFact;
    for n = gene_freq_histogram_observed(1)+1 : N
        fact = fact + log(n);
    end
    
    % optimize gclass mixture coefficients and probs assuming N genes in supragenome
    params = fmincon( @likelihood, ML_params, A, b, Aeq, beq, lb, ub, [], options );
    
    % normalize log probability
    logProb = -likelihood( params ) + fact;
    
    % add result to X array
    X = [X; [N logProb] ];
    
    % check if this is the max likelihood
    if ( logProb > ML_logProb )
        ML_logProb = logProb;
        ML_params = params;
        ML_N = N;
    end
    
    N+1;
end

dlmwrite(fileName_N_vs_likelihood, X, 'delimiter', '\t');

%Echo results (notice the absence of semi-colon line endings):
ML_N
ML_logProb
ML_gclass_mixcoeff = ML_params(1:K)
ML_gclass_probs = ML_params(K+1:end)

run_stats(...
    ML_N,...
    ML_gclass_mixcoeff,...
    ML_gclass_probs,...
    gene_freq_histogram_predicted_strains,...
    new_core_genes_predicted_strains...
    );

%{
**************************************************************************
SECTION 4.  FUNCTIONS:
            get_prompts()
            get_run_option()
            get_input_numeric(default_numeric, parameter_description)
            likelihood(params)
            run_stats(...
                ML_N,...
                ML_gclass_mixcoeff,...
                ML_gclass_probs,...
                gene_freq_histogram_predicted_strains,...
                new_core_genes_predicted_strains...
                )

**************************************************************************
    %}
    function get_prompts()
        prompt_fileName_diary = [...
            'Enter a string (e.g., ''Saureus17'', without quotes) to yield a ',...
            'run-specific name (''CommandWindow_Saureus17.txt'' in this example) ',...
            'for your MatLab session''s diary file.\n'...
            ];
        prompt_fileName_N_vs_likelihood = [...
            'Enter a string (e.g., ''Saureus17'', without quotes) to yield a ',...
            'run-specific name (''N_vs_Likelihood_Saureus17.txt'' in this example) ',...
            'for your ''N_vs_likelihood'' output file.\n'...
            ];
        prompt_vector_entry = [...
            'Vectors are entered inside square parentheses with ',...
            'their elements separated by spaces , e.g.: [e1 e2 e3 e4]\n'
            ];
        prompt_gene_freq_histogram_observed = [...
            'Enter the gene frequency histogram vector for the observed ',...
            'number of genes in n = 0, 1, 2, ...,|S|-1, |S| strains examined. ',...
            prompt_vector_entry
            ];
        prompt_gclass_probs_start = [...
            '(K) initial values of the (mu) vector for the (population) gene class',...
            'frequenies.  This vector''s final ML values will be estimated during ',...
            'the optimization. The K-th element of this vector should be equal to ',...
            '1.00 (for the core genes). ',...
            prompt_vector_entry...
            ];
        prompt_gclass_probs_separation_min = [...
            'value of the minimum separation between gene-class frequencies used ',...
            'during the optimization.\n'...
            ];
        prompt_gclass_probs_min = [...
            'value of the minimum gene-class frequency used ',...
            'during the optimization.\n'...
            ];
        prompt_gclass_mixcoeff_min = [...
            'value of the minimum mixture coefficient used ',...
            'during the optimization.\n'...
            ];
        prompt_N_start_offset = [...
            'value of the ''offset'' (to add to the total number of observed ',...
            'genes) to yield ''Nmin'', the start of the range (from Nmin to Nmax) ',...
            'of supragenome size (N) values to examine during the optimization.\n'...
            ];
        prompt_N_range = [...
            'value (Nmax - Nmin) of the range of supragenome size (N) values to ',...
            'examine during the optimization.\n'...
            ];
        prompt_gene_freq_histogram_predicted_strains = [...
            'Enter a value for the number of strains |S| in the predicted gene ',...
            'frequency histogram vector (the predicted number of genes in ',...
            'n = 0, 1, 2, ...,|S|-1, |S| strains).  These predictions will be ',...
            'made using the ML estimates of N (supragenome size), mu (gene class ',...
            'population frequencies), and pi (gene class mixture coefficients).\n'...
            ];
        prompt_new_core_genes_predicted_strains = [...
            'value for the number of strains for which new and core genes will ',...
            'be predicted.  These predictions will be made using the ML ',...
            'estimates of N (supragenome size), mu (gene class population ',...
            'frequencies), and pi (gene class mixture coefficients).\n'...
            ];
        prompt_ML_N = [...
            'Enter the maximum likelihood value for N (supragenome size).\n'...
            ];
        prompt_ML_gclass_mixcoeff = [...
            'Enter the maximum likelihood values for the pi vector of gene class ',...
            'mixture coefficients. ',...
            prompt_vector_entry
            ];
        prompt_ML_gclass_probs = [...
            'Enter the maximum likelihood values for the mu vector of gene class ',...
            'probabilities (population frequencies). ',...
            prompt_vector_entry
            ];
    end

    function run_option = get_run_option()
        prompt_run_option = [...
            'Enter ''1'' or ''2'' (without quotes) to select a run option:\n',...
            'NOTE: You may also select the default run option [1] with a carriage return.\n\n',...
            '[1] Run the optimization and the statistics.\n',...
            ' 2  Run the statistics only.\n'...
            ];
        run_option_buffer = input(prompt_run_option);
        if isempty(run_option_buffer)
            run_option = 1;
        else
            run_option = int32(run_option_buffer);
        end
    end

    function input_numeric = get_input_numeric(default_numeric, parameter_description)
        input_prompt = [...
            'Press return to accept the default shown in parentheses, or enter an ',...
            'alternative for the ',...
            parameter_description,...
            '[', num2str(default_numeric), ']'...
            ];
        input_buffer = input(input_prompt);
        %If 'input' receives a carriage return, the result is an empty matrix
        if isempty(input_buffer)
            input_numeric = default_numeric;
        else
            input_numeric = input_buffer;
        end
    end

    function like = likelihood(params)
        condP = zeros(S+1,K);
        for i = 0:S
            for k=1:K
                condP(i+1,k) = nchoosek(S,i) * params(K+k)^i * ( 1-params(K+k) )^(S-i);
            end
        end
        like = 0;
        for i = 0:S
            temp = 0;
            for k = 1:K
                temp = temp + params(k)*condP(i+1,k);
            end
            like = like + gene_freq_histogram_observed(i+1) * log(temp);
        end
        % we need to maximize using the minimize function, so take negative.
        like = -real(like);
    end

    function run_stats(...
            ML_N,...
            ML_gclass_mixcoeff,...
            ML_gclass_probs,...
            gene_freq_histogram_predicted_strains,...
            new_core_genes_predicted_strains...
            )
        K = length(ML_gclass_mixcoeff);
        S = gene_freq_histogram_predicted_strains;
        
        % generate histogram of the predicted number of genes in n of
        % |S| (= gene_freq_histogram_predicted_strains) strains examined.
        gene_freq_histogram_predicted = zeros(S,1);
        for n=1:S
            for k=1:K
                gene_freq_histogram_predicted(n) = gene_freq_histogram_predicted(n) + ...
                    nchoosek(S,n) * ML_gclass_probs(k)^n * (1-ML_gclass_probs(k))^(S-n)...
                    * ML_gclass_mixcoeff(k);
            end
        end
        
        % generate predictions of core genes and new genes
        core = zeros(new_core_genes_predicted_strains,1);
        core_var = zeros(new_core_genes_predicted_strains,1);
        new = zeros(new_core_genes_predicted_strains,1);
        new_var = zeros(new_core_genes_predicted_strains,1);
        total = sum(ML_gclass_mixcoeff(1:end))*ones(new_core_genes_predicted_strains,1);
        total_var = zeros(new_core_genes_predicted_strains,1);
        
        for n=1:new_core_genes_predicted_strains
            for k = 1:K
                core(n) = core(n) + ML_gclass_probs(k)^n * ML_gclass_mixcoeff(k);
                core_var(n) = core_var(n) + ML_gclass_probs(k)^n * ...
                    (1-ML_gclass_probs(k)^n) * ML_gclass_mixcoeff(k);
                new(n) = new(n) + ML_gclass_probs(k)*...
                    (1-ML_gclass_probs(k))^(n-1)*ML_gclass_mixcoeff(k);
                new_var(n) = new_var(n) + ...
                    ML_gclass_probs(k)*(1-ML_gclass_probs(k))^(n-1) * ...
                    (1-ML_gclass_probs(k)*(1-ML_gclass_probs(k))^(n-1)) * ...
                    ML_gclass_mixcoeff(k);
                total(n) = total(n) - (1-ML_gclass_probs(k))^n * ML_gclass_mixcoeff(k);
                total_var(n) = total_var(n) + (1-ML_gclass_probs(k)^n) * ...
                    ML_gclass_probs(k)^n * ML_gclass_mixcoeff(k);
            end
        end
        
        % calculate pairwise shared and different genes
        shared = 0;
        shared_var = 0;
        different = 0;
        different_var = 0;
        
        for k = 1:K
            shared = shared + ML_gclass_probs(k)^2 * ML_gclass_mixcoeff(k);
            shared_var = shared_var + ML_gclass_probs(k)^2 * ...
                (1 - ML_gclass_probs(k)^2) * ML_gclass_mixcoeff(k);
            different = different + 2*ML_gclass_probs(k)*(1-ML_gclass_probs(k)) * ...
                ML_gclass_mixcoeff(k);
            different_var = different_var + ...
                2*ML_gclass_probs(k)*(1-ML_gclass_probs(k)) * ...
                (1 - 2*ML_gclass_probs(k)*(1-ML_gclass_probs(k))) * ...
                ML_gclass_mixcoeff(k);
        end
        
        % calculate total genes per genome
        strain_total = 0;
        strain_total_var = 0;
        for k = 1:K
            strain_total = strain_total + ML_gclass_probs(k) * ML_gclass_mixcoeff(k);
            strain_total_var = strain_total_var + ML_gclass_probs(k) * ...
                (1-ML_gclass_probs(k)) * ML_gclass_mixcoeff(k);
        end
        
        %Echo results (notice the absence of semi-colon line endings):
        gene_freq_histogram_predicted = int32(gene_freq_histogram_predicted * ML_N)
        
        core = int32(core * ML_N)
        new = int32(new * ML_N)
        total = int32(total * ML_N)
        
        core_stdv = int32(((core_var).^(1/2)) * ML_N)
        new_stdv = int32(((new_var).^(1/2)) * ML_N)
        total_stdv = int32(((total_var).^(1/2)) * ML_N)
        
        shared = int32(shared * ML_N)
        shared_stdv = int32(((shared_var)^(1/2)) * ML_N)
        
        different = int32(different * ML_N)
        different_stdv = int32(((different_var)^(1/2)) * ML_N)
        
        strain_total = int32(strain_total * ML_N)
        strain_total_stdv = int32((strain_total_var^(1/2)) * ML_N)
        
        %Echo execution data (notice the absence of semi-colon line endings):
        programStart
        programFinished = datestr(now)
    end
end

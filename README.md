# Kappa rules for gene transcription factors binding 
   Python code based on PySb (https://github.com/pysb/pysb) to generate Kappa rules for n-sites Gene-Transcription factor binding combination

# Input: 
        python GeneTranscriptionFactorKappaRule -i #transcription_factors -o <output_file>

# Output:

- For #transcription_factors=2


-------------------------------------- Kappa rules ------------------------------------


% agent: TF(a_0,a_1)

% agent: G(t_0,t_1)


%var: 'k_fwd_1' 1.000000e-02

%var: 'k_fwd_2' 1.000000e-02

%var: 'k_fwd_3' 1.000000e-02

%var: 'k_fwd_4' 1.000000e-02

'binding_rule_1' G(t_0,t_1),TF(a_0,a_1) -> G(t_0!0,t_1),TF(a_0!0,a_1) @ 'k_fwd_1'

'binding_rule_2' G(t_0,t_1),TF(a_0,a_1) -> G(t_0,t_1!1),TF(a_0,a_1!1) @ 'k_fwd_2'

'binding_rule_3' G(t_0,t_1!1),TF(a_0,a_1!1) -> G(t_0!0,t_1!1),TF(a_0!0,a_1!1) @ 'k_fwd_3'

'binding_rule_4' G(t_0!0,t_1),TF(a_0!0,a_1) -> G(t_0!0,t_1!1),TF(a_0!0,a_1!1) @ 'k_fwd_4'

-------------------------------------- ########### ------------------------------------


# Pictorial overview:

                                    	  _  _  (Gene)
					 	
                                   		 /    \
                                   		/      \
 			 
 		(binding_rule_1)	     *_ _      _ *_        (binding_rule_2)
				
                                       	\       /
                                           	 \     / 
						 
                        	(binding_rule_4)  *_ _*  (binding_rule_3)


                                


# GeneTranscriptionFactorBindingKappaRule
Python code to generate Kappa rules for an n-site Gene-Transcription factor binding combination

#Input: agent=['Gene','TranscriptionFactor']
        site_type={'G':['t_'], 'TF':['a_']}
        _number_of_tf=n (number of transcripton factors)
		number_of_sites ={'G':[_number_of_tf], 'TF':[_number_of_tf]} 

#Output:

For n=2

-------------------------------------- Kappa rules ------------------------------------
%agent: TF(a_0,a_1)
%agent: G(t_0,t_1)

%var: 'k_fwd_1' 1.000000e-02
%var: 'k_fwd_2' 1.000000e-02
%var: 'k_fwd_3' 1.000000e-02
%var: 'k_fwd_4' 1.000000e-02

'binding_rule_1' G(t_0,t_1),TF(a_0,a_1) -> G(t_0!0,t_1),TF(a_0!0,a_1) @ 'k_fwd_1'
'binding_rule_2' G(t_0,t_1),TF(a_0,a_1) -> G(t_0,t_1!1),TF(a_0,a_1!1) @ 'k_fwd_2'
'binding_rule_3' G(t_0,t_1!1),TF(a_0,a_1!1) -> G(t_0!0,t_1!1),TF(a_0!0,a_1!1) @ 'k_fwd_3'
'binding_rule_4' G(t_0!0,t_1),TF(a_0!0,a_1) -> G(t_0!0,t_1!1),TF(a_0!0,a_1!1) @ 'k_fwd_4'
-------------------------------------- ########### ------------------------------------

#Pictorial overview:

                                    	 _ 	_  (Gene)
                                     	t_0 t_1
                                     
                                   		  /    \
                                   		 /      \
 (binding_rule_1)			   (a_0)  *_ _      _ *_ (a_1)     		(binding_rule_2)
                                         \       /
                                          \     / 
                         (binding_rule_4)  *_ _*  (binding_rule_3)


                                



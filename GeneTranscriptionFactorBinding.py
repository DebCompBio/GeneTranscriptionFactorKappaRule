#/usr/bin/python
# -*- coding: utf-8 -*-
"""
:mod:`tf_binding_kappa` -- Kappa rules for transcription factor binding

                            - - - - -
                            /   |   \  .........
                           /    |    \
                     _ *_ _ _ _  _ _ *_ _ _  ......
                       /  \            \
                      /    \            \
            _ *_ _ *_ _   _ *_ _ _ *_   .........
                     \                       /
                      \         :::         /
                                           /
                       :      
                           *_ *_ *_ *_ *_ Pol
                                 
                        Transcription initiation complex
* : TF factors
- : TF binding sites 
===================================
.. module:: tf_binding_kappa
 :platform: Unix, Windows
 :synopsis: Generate .ka script of rules for n transcrption factors binding.
.. moduleauthor:: Debdas Paul
Requirement(s)::
 1. You have to set the KaSim excutable path by using pysb.kappa.set_kappa_path('KaSim path') 
"""

import pysb
from pysb import *
from pysb.export import export
from pysb import kappa
import copy
import operator
import numpy as np
import re
from operator import itemgetter
import itertools
""" Setting the KaSim path 

""" 
pysb.kappa.set_kappa_path('./KaSim')

#Global
Model()


class GeneTranscriptionFactorBinding(object):

    def __init__(self,agent,site_type,number_of_sites):

      self.agent = agent
      self.site_type = site_type
      self.number_of_sites = number_of_sites

    def _site_enumeration(self):
      """  

          Return: 

      """
      site_type_dict = {}  # site_type_dict = {'TF/G_sites_type_x_': ['x_k']}

      for key in site_type.keys():
          for length in range(len(site_type[key])): 
             site_type_dict[str(key)+'_sites_type_'+str(site_type[key][length])] = [site_type[key][length] + '%d' % i for i in range(number_of_sites[key][length])] 
      
      return site_type_dict 

    def _creating_monomers(self, e_s,s_t):
        """  
          
          Arg: 


          Return: 

      """
        for site_t in s_t.keys():
          concat_sites = []
           
          for site_e in e_s.keys():
            if re.split('\_', site_e)[0]==str(site_t):
               concat_sites=concat_sites + e_s[site_e]

                 
          Monomer(site_t, concat_sites)


    def _bounded_unbounded(self,e_s,s_t):  
        _unbounded_sites ={}
        _bounded_sites   = {}
        for e_s_key in e_s.keys(): 
         _un_bound_key=re.split('\_', e_s_key)[0]+'_'+re.split('\_', e_s_key)[3]+'_'  

         _unbounded_sites[_un_bound_key]= dict([(site_name, None) for site_name in e_s[e_s_key]])
             
         _bounded_sites[_un_bound_key] = dict([(site_name, bond_num)
               for bond_num, site_name in enumerate(e_s[e_s_key])])
        

        return (_unbounded_sites, _bounded_sites)

    def _powerset_generation(self, tf):
   
        result = [[]]
        for x in sorted(tf):
             result.extend([subset + [x] for subset in result])
        result = result[1:]
        return sorted(result, key=len)
    def _sorting(self,ss):
        return sorted(ss, key=len)

    def _product_level(self,sorted_dic,_b_sites,_b_sites_original,agent_name):
        
        _pro_level={}
        for dic_iter in range(1,_number_of_tf+1):
             _pro_level[dic_iter]=[]

        for _pro_key in _pro_level.keys():
            for sorted_dic_length in range(len(sorted_dic)):
                if int(_pro_key) == len(sorted_dic[sorted_dic_length]):
                   
                   for _b_site_key in _b_sites.keys():
                       if  re.split('\_',_b_site_key)[0]==agent_name: 
                          for _b_site_key_key in _b_sites[_b_site_key].keys():
                          
                                if _b_site_key_key  not in  sorted_dic[sorted_dic_length]:
                                   _b_sites[_b_site_key][_b_site_key_key]=None
                                     
                          _pro_level[_pro_key].append(_b_sites[_b_site_key]) 
                          _b_sites=copy.deepcopy(_b_sites_original)
        return  _pro_level

    def _reverse_product_level(self,_pro_level):
         
        _pro_level_reverse=copy.deepcopy(_pro_level) 

        for _rev_pro_key in _pro_level_reverse.keys():
            for _rev_pro_key_length in range(len(_pro_level_reverse[_rev_pro_key])):
                for _rev_pro_key_key in _pro_level_reverse[_rev_pro_key][_rev_pro_key_length].keys():
                    if _pro_level_reverse[_rev_pro_key][_rev_pro_key_length][_rev_pro_key_key]==None:
                       _pro_level_reverse[_rev_pro_key][_rev_pro_key_length][_rev_pro_key_key]= int(re.split('\_',_rev_pro_key_key)[1])
                    else:
                       _pro_level_reverse[_rev_pro_key][_rev_pro_key_length][_rev_pro_key_key] = None
                  
        return _pro_level_reverse
    def _generate_kappa_rules(self,rule_index,rule_iter,
                                _pro_level_tf,_rev_pro_level_tf,_pro_level_gene,_rev_pro_level_gene):
        """First level rules

        """ 
        for iter_n in range(_number_of_tf,_number_of_tf-1,-1):
               for iter_pro_level in range(0,len(_rev_pro_level_gene[iter_n])): 
                  for iter_rev_pro_level in range(0,len(_pro_level_gene[rule_iter])):  
                    Rule ('binding_rule_'+str(rule_index),
                      G(**_rev_pro_level_gene[iter_n][iter_pro_level]) + TF(**_rev_pro_level_tf[iter_n][iter_pro_level]) >>
                      G(**_pro_level_gene[rule_iter][iter_rev_pro_level]) + TF(**_pro_level_tf[rule_iter][iter_rev_pro_level]),Parameter('k_fwd_'+str(rule_index), 1e-2))
                    rule_index=rule_index+1 
               rule_iter=rule_iter+1 
        if _number_of_tf >1:

            """ Intermediate rules
            
            """              
            true_flag=0 
            for iter_n in range(_number_of_tf-1,1,-1):
                for iter_pro_level in range(0,len(_rev_pro_level_gene[iter_n])): 
                    for iter_rev_pro_level in range(0,len(_pro_level_gene[rule_iter])): 
                        
                        for key in _rev_pro_level_gene[iter_n][iter_pro_level].keys():
                            if _pro_level_gene[rule_iter][iter_rev_pro_level][key] == _rev_pro_level_gene[iter_n][iter_pro_level][key]:
                                true_flag=true_flag+1 

                        if true_flag==_number_of_tf-1:
                           Rule ('binding_rule_'+str(rule_index),
                                G(**_rev_pro_level_gene[iter_n][iter_pro_level]) + TF(**_rev_pro_level_tf[iter_n][iter_pro_level]) >>
                                G(**_pro_level_gene[rule_iter][iter_rev_pro_level]) + TF(**_pro_level_tf[rule_iter][iter_rev_pro_level]),Parameter('k_fwd_'+str(rule_index), 1e-2))
                           rule_index=rule_index+1
                        true_flag=0 
                rule_iter=rule_iter+1 

            """ Last level rules
            
            """    
            for iter_n in range(1,0,-1):
               for iter_pro_level in range(0,len(_rev_pro_level_gene[iter_n])): 
                  for iter_rev_pro_level in range(0,len(_pro_level_gene[rule_iter])):  
                    Rule ('binding_rule_'+str(rule_index),
                      G(**_rev_pro_level_gene[iter_n][iter_pro_level]) + TF(**_rev_pro_level_tf[iter_n][iter_pro_level]) >>
                      G(**_pro_level_gene[rule_iter][iter_rev_pro_level]) + TF(**_pro_level_tf[rule_iter][iter_rev_pro_level]),Parameter('k_fwd_'+str(rule_index), 1e-2))
                    rule_index=rule_index+1 
               rule_iter=rule_iter+1    
_number_of_tf=1         
agent=['Gene','TranscriptionFactor']
site_type={'G':['t_'], 'TF':['a_']}
number_of_sites ={'G':[_number_of_tf], 'TF':[_number_of_tf]}
GtObj = GeneTranscriptionFactorBinding(agent, site_type, number_of_sites)
enumerate_site=GtObj._site_enumeration()
GtObj._creating_monomers(enumerate_site, site_type)
[_ub_sites, _b_sites]=GtObj._bounded_unbounded(enumerate_site,site_type)
transcription_factor_power_set=GtObj._powerset_generation(_ub_sites['TF_a_'].keys())
gene_power_set=GtObj._powerset_generation(_ub_sites['G_t_'].keys())
_b_sites_original=copy.deepcopy(_b_sites)
_pro_level_tf=GtObj._product_level(transcription_factor_power_set,_b_sites,_b_sites_original,'TF')
_b_sites=copy.deepcopy(_b_sites_original)
_pro_level_gene=GtObj._product_level(gene_power_set,_b_sites,_b_sites_original,'G')
_rev_pro_level_tf=GtObj._reverse_product_level(_pro_level_tf)
_rev_pro_level_gene=GtObj._reverse_product_level(_pro_level_gene)
rule_index    =1
rule_iter     =1
GtObj._generate_kappa_rules(rule_index,rule_iter,_pro_level_tf,_rev_pro_level_tf,_pro_level_gene,_rev_pro_level_gene)
kappa_str = export(model, 'kappa')
print(kappa_str)

""" Writing the rules to a file
"""
#with open('GeneTranscriptionFactorBinding.ka', 'wt') as f:
#       f.write(kappa_str)








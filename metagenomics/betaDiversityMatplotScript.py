# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 16:57:29 2019

@author: derek.bickhart-adm
"""

import argparse as ap
import math
import sys, os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from sklearn import decomposition
from sklearn import preprocessing as prep
from sklearn import manifold
import matplotlib.gridspec as gridspec
from sklearn.metrics import pairwise_distances
import seaborn as sns
import pandas as pd
import numpy as np
import itertools
from scipy import stats

def arg_parse():
    parser = ap.ArgumentParser(
            description = "A rewrite of the cumbersome betadiversity script"
            )
    parser.add_argument('-t', '--table', 
                        help="Input tab-delimited OTU count table with sample classification as the last column",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="output base name",
                        required=True, type=str,
                        )
    parser.add_argument('-f', '--featureID',
                        help="The name of column with the feature ID",
                        type=str, default="class"
                        )
    return parser.parse_args()

class betaData(object):
    
    def __init__(self, args, file : str, featureId : str):
        # one-off functions
        self.is_a_comb = lambda c, _l_ : bool(c in _l_) 
        self.they_overlap = lambda pos_a1, pos_a2, pos_b1, pos_b2 : bool( min([pos_a2, pos_b2]) - max([pos_a1, pos_b1]) )
        self.p_value = lambda group_a, group_b : stats.ttest_ind(group_a, group_b, axis=0, equal_var=False)[1]
        
        # Hard coded settings that I can change later
        self.args = vars(args)
        self.args['algorithm'] = "mds"
        self.args['distance'] = 'braycurtis'
        self.shapes_types = '.,ov^<>1234sp*hH+xDd|_'
        self.args['facecolor'] = 'white'
        self.args['p_values_only_for'] = []
        self.args['boxplot'] = True
        
        self.raw = pd.read_csv(file, sep='\t', header=0)
        
        # separate the data from the features
        self.data = self.raw.drop(featureId, axis=1)
        self.features = self.raw[featureId].tolist()
        
    def betaDiversity(self):
        """
        Main calculation routine for betadiversity
        Returns:
           transformed data
           exp_var
        """
        return self.transformation(self.args['algorithm'], self.data, self.args['distance'])

    def transformation(self, way, data, distance):
        T_func = (self.pca if way == 'pca' else (self.mds if way == 'mds' else self.nmds)) if not self.args['boxplot'] else 'just_the_matrix'
        if T_func != 'just_the_matrix':
            if T_func == self.pca: 
                transformed, exp_var = T_func(data)
            else: 
                transformed, exp_var = T_func(self.compute_distances(data,distance)), None
            return transformed, exp_var
        else:
            return self.compute_distances(data,distance), None        
        
    def pca(self, f):
        pca = decomposition.PCA(n_components=2)
        return pca.fit_transform(f), pca.explained_variance_ratio_

    def mds(self, d):
        try:
            mds = manifold.MDS(n_components=2, max_iter=5000, eps=1e-9, dissimilarity='precomputed')
            return mds.fit(d).embedding_
        except ValueError:
            print('You Have NaNs in the data: here the coordinates')
            exit(1)     

    def nmds(self, d):
        nmds = manifold.MDS(n_components=2, metric=False, max_iter=3000, eps=1e-12, dissimilarity="precomputed", n_jobs=1, n_init=1)
        return nmds.fit_transform(d)

    def compute_distances(self, data, metric):
        if metric == 'precomputed': return data
        elif metric == 'lbraycurtis':
            ldata = np.matrix([[(math.log(1.0+float(l)) if float(l) > 0.0 else 0.0) for l in d] for d in data])
            return pairwise_distances(ldata, metric='braycurtis')
        elif metric == 'sbraycurtis':
            sdata = np.matrix([[(math.sqrt(float(l)) if float(l) > 0.0 else 0.0) for l in d] for d in data])
            return pairwise_distances(sdata, metric='braycurtis')
        else:
            return pairwise_distances(data, metric='braycurtis')     
        
        def save_coordinates(self, coordinates_f): self.sample_and_coordinates.to_csv(coordinates_f, sep='\t', header=True, index=True)
        
    def betadiv_statistical_support(self, text):
        print('\n=============================================')
        print('Statistical Support ' + ('about ' + self.args['title'] if self.args['title'] else ''))     
        print('=============================================\n')
        for group_a,group_b in self.combinations_of_classes:
            statistic,p_value = stats.ttest_ind(self.groups[group_a], self.groups[group_b], axis=0, equal_var=False)
            print(' - '.join(list((group_a,group_b)))) 
            print(' P-PVALUE = ' + str(p_value))
            print('sign. TRUE' if (p_value < 0.05) else '')
        print('\n=============================================\n')
        
    def p_value_on_plot(self, cup, drop, level, just_zerotwo=False):
        if not self.done_[cup]:  
            self.done_[cup] = True
            p_ = self.p_value(self.groups[cup[0]], self.groups[cup[1]])    
            if p_ < 0.05:
                plt.plot(\
                   [self.indices[cup[0]] + ((drop if not just_zerotwo else just_zerotwo) if self.indices[cup[0]] < self.indices[cup[1]] else -(drop if not just_zerotwo else just_zerotwo))\
                  , self.indices[cup[0]] + ((drop if not just_zerotwo else just_zerotwo) if self.indices[cup[0]] < self.indices[cup[1]] else -(drop if not just_zerotwo else just_zerotwo))\
                  , self.indices[cup[1]] - ((drop if not just_zerotwo else just_zerotwo) if self.indices[cup[0]] < self.indices[cup[1]] else -(drop if not just_zerotwo else just_zerotwo))\
                  , self.indices[cup[1]] - ((drop if not just_zerotwo else just_zerotwo) if self.indices[cup[0]] < self.indices[cup[1]] else -(drop if not just_zerotwo else just_zerotwo))]\
                 , [level, level + drop*0.5, level + drop*0.5, level] if self.args['p_values'] == 'above' \
                    else [level, level - drop*0.5, level - drop*0.5, level], lw=1., color='k')
                
                plt.text((self.indices[cup[0]]+self.indices[cup[1]])*0.5, level\
                        ,'p. '+str(p_), ha='center', va='bottom', color='k', fontsize=3)
                return p_
        return False
    
    def box_plot(self, dist_mat):       
        sns.set_style(self.args['facecolor'])  
        #for c in self.couples_by_class: self.couples_by_class[c], 'Apparently, I cover one thousand whales!'

        print(self.args['p_values_only_for'])

        class box_plot_object(object):
            def __init__(self, class_, color_, cat_var_, couple_of_samples, dist_mat):
                self.class_ = class_
                self.color_ = color_
                self.cat_var_ = cat_var_
                # The couple_of_samples[0] is the class and couple_of_samples[1] is the sample index
                self.beta_diversity = dist_mat.get_value(couple_of_samples[0], couple_of_samples[1])
                
        tcolors = ["crimson", "blue", "green"]
        data = pd.DataFrame(\
                [[ob.class_, ob.color_, ob.cat_var_, ob.beta_diversity]   \
                  for ob in [\
                  box_plot_object(cl,co,ct,cp,self.dist_mat) for cl,co,ct,cp in zip(\
                    list(itertools.chain.from_iterable(\
                       [[c[0] for i in range(len(self.couples_by_class[c[0]]))] for c in self.legend if len(self.couples_by_class[c[0]]) ]))\
                  , list(itertools.chain.from_iterable(\
                       [[c[1] for i in range(len(self.couples_by_class[c[0]]))] for c in self.legend if len(self.couples_by_class[c[0]]) ]))\
                  , list(itertools.chain.from_iterable(\
                       [[c[2] for i in range(len(self.couples_by_class[c[0]]))] for c in self.legend if len(self.couples_by_class[c[0]]) ]))\
                  , list(itertools.chain.from_iterable(\
                       [[couple for couple in self.couples_by_class[c[0]]] for c in self.legend if len(self.couples_by_class[c[0]])]))\
                  ) ]], columns=['', 'color', 'group_by', 'Beta-Diversity']) 


        ax = sns.swarmplot(data=data, x='', y='Beta-Diversity'\
          , hue=None if len(list(set(data['group_by'].tolist())))==1 else 'group_by', dodge=True, s=2, color='black', alpha=1.)
        ax = sns.boxplot(data=data, x='', y='Beta-Diversity'	
          , hue=None if len(list(set(data['group_by'].tolist())))==1 else 'group_by', palette=dict([(c[0], c[1]) for c in self.legend]))

        data.columns = ['class', 'color', 'group_by', 'Beta-Diversity']
        self.groups = dict([(cl, data.loc[data['class'].isin([cl]), 'Beta-Diversity'].tolist()) for cl in list(set(data['class'].tolist()))])
        self.combinations_of_classes = list(map(tuple, itertools.combinations(list(set(data['class'].tolist())), 2)))
               
        self.distributions = []
        for class_ in data['class'].tolist(): 
            if class_ not in self.distributions: self.distributions.append(class_) 

        self.indices = dict([(d,e) for e,d in enumerate(self.distributions)])
        self.couple_levels = dict()
        drop = data['Beta-Diversity'].max()*0.02

        if self.args['p_values'] == 'above': H = data['Beta-Diversity'].max() + drop #* 2
        else: H = data['Beta-Diversity'].min() - drop ##* 2

        level = H
        #cover = dict([(dat,level) for dat in self.distributions])
        distances = dict([(couple, (max([self.indices[couple[0]], self.indices[couple[1]]]) \
                                  - min([self.indices[couple[0]], self.indices[couple[1]]]))) \
                                    for couple in self.combinations_of_classes])
        self.done_ = dict([(couple, False) for couple in self.combinations_of_classes])
        ##prevs = {} #min(distances.values())

        distance = 1
        for i in range(1, len(self.distributions)): ##list(itertools.chain.from_iterable(range(prev, len(self.distributions)), ):
            cup = self.distributions[i], self.distributions[i-1]   
            if self.is_a_comb(cup, self.combinations_of_classes): print('coupling OK')
            else:
                cup = (cup[1],cup[0])
                if self.is_a_comb(cup, self.combinations_of_classes): print('reversed coupling OK')
            sign = self.p_value_on_plot(cup, drop, level)
        
        level = (level + drop * 2) if self.args['p_values'] == 'above' else (level - drop * 2)
        #print distances, '  guarda questo'

        for distance in range(2, max(distances.values())+1):
            for any_comb in self.combinations_of_classes:
                if distances[any_comb] == distance: 
                   sign =  self.p_value_on_plot(any_comb, drop, level)  
                   if sign: 
                       level = (level + drop * 2) if self.args['p_values'] == 'above' else (level - drop * 2 )   

        plt.setp(ax.get_xticklabels(), rotation=38, ha='right')
        plt.subplots_adjust(bottom=.3) 
        plt.suptitle(self.args['title'], fontsize=10)
        plt.savefig(self.args['stdout']+'.'+self.args['format'], dpi=400)

        self.betadiv_statistical_support(self.args['text_on'])

def main(args):
    args = arg_parse()
    bd = betaData(args, args.table, args.featureID)
    
    transformed_data, expr_var = bd.betaDiversity()
    
    transformed_data
    expr_var
    
if __name__ == "__main__":
    args = arg_parse()
    main(args)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 15:08:18 2023

@author: jwahnzavalet
"""
import json
import requests
import numpy as np

class DGIAPI:
    'API Example class for DGI API.'
    domain = 'https://dgidb.org/'
    api_path = '/api/v1/interactions.json'
    
    def __init__(self, genes=[], interaction_sources=[], interaction_types=[], 
                 gene_categories=[], source_trust_levels=[], antineoplastic_only=[],
                 save=True):

        self.genes = genes
        self.interaction_sources = interaction_sources
        self.interaction_types = interaction_types
        self.gene_categories = gene_categories
        self.source_trust_levels = source_trust_levels
        self.antineoplastic_only = antineoplastic_only
        self.save = save
        
    def run_workflow(self):
        self.create_request()
        self.post_request()
        if (not self.save):
            self.print_response()
        else:
            return(self.compile_df())
        
    def create_request(self):
        # Initialize the payload dictionary
        self.payload = {}
        
        # Define a dictionary of variables and their corresponding attributes
        variables = {
            'genes': self.genes,
            'interaction_sources': self.interaction_sources,
            'gene_categories': self.gene_categories,
            'interaction_types': self.interaction_types,
            'source_trust_levels': self.source_trust_levels
        }
        
        # Iterate through the variables and their values
        for var, value in variables.items():
            # Check if the value is a list before applying join
            if isinstance(value, list):
                value = [item.replace(".", " ") for item in value]
                self.payload[var] = ",".join(value)
            else:
                if (value != 'None') & (value is not None):
                    value = value.replace(".", " ")
                self.payload[var] = value
            
        # Add antineoplastic drug type if applicable
        if self.antineoplastic_only is not None:
            self.payload['drug_types'] = 'antineoplastic'
        
        print(self.payload)
        
        
    def post_request(self):
        self.request = DGIAPI.domain + DGIAPI.api_path
        self.response = requests.post(self.request, data = self.payload)
    def compile_df(self):
        data = []
        response = json.loads(self.response.content.decode('utf-8'))
        
        matches = response['matchedTerms']
        
        if matches:
            for match in matches:
                gene = match['geneName']
                categories = match['geneCategories']
                categories.sort()
                joined_categories = ",".join(categories)
                for interaction in match['interactions']:
                    source = interaction['source']
                    drug = interaction['drugName']
                    interaction_type = interaction['interactionType']
                    data.append({
                        'gene_name': gene,
                        'drug_name': drug,
                        'interaction_type': interaction_type,
                        'source': source,
                        'gene_categories': joined_categories.lower()
                    })
        if data:
            df = np.matrix(data)
            return df
        else:
            return 
        
    def print_response(self):
        response = json.loads(self.response.content.decode('utf-8'))
        matches = response['matchedTerms']
        if(matches):
          print ("gene_name\tdrug_name\tinteraction_type\tsource\tgene_categories")
        for match in matches:
            gene = match['geneName']
            categories = match['geneCategories']
            categories.sort()
            joined_categories = ",".join(categories)
            for interaction in match['interactions']:
                source = interaction['source']
                drug = interaction['drugName']
                interaction_type = interaction['interactionType']
                print (gene + "\t" + drug + "\t" + interaction_type + "\t" + source + "\t" + joined_categories.lower())
                
        for unmatched in response['unmatchedTerms']:
            print ("Unmatched search term: " + unmatched['searchTerm'])
            print ("Possible suggestions: " + ",".join(unmatched['suggestions']))

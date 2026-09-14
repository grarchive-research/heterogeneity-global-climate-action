"""
@author: Giorgio Ricchiuti
@website: www.grarchive.net

Description
-----------
Loads country-level weights and emissions used by the agent-based
model (Section 4) to introduce heterogeneity in country size and
emissions: countries' relative GDP share (weights) and their GHG
emissions, as country-name-keyed arrays aligned by country.

Data source: countries_size_emissions.csv (World Bank GDP shares,
1989-2019 average; see main paper for details)

Version 3.0, Jan-Feb 2025
"""

import numpy as np
import pandas as pd

def dataset(): 
    
    file_path = "countries_size_emissions.csv"  
    data = pd.read_csv(file_path)
    #print(data.head())
    
    country_code = data["countrycode"]  # Get the list of countries
    country_names = data["CountryName"]  # Assuming column "CountryName" exists
    
    weights = data["weights"].dropna()  # Remove NaN values
    GHG=data["GHG"].dropna()
    
    min_weight = weights.min()
    max_weight = weights.max() # Get the weight values
    weight_dict = data.set_index("countrycode")["weights"].to_dict()
    
    min_GHG = GHG.min()
    max_GHG = GHG.max() # Get the weight values
    GHG_dict = data.set_index("countrycode")["GHG"].to_dict()
    
    # Build country-name-keyed dicts, then convert to arrays in the same
    # order as country_names (used by the ABM to weight countries by
    # size/emissions rather than treating them as identical)
    weights_array = np.array(weights)
    country_weight_dict = {country_name: weight for country_name, weight in zip(country_names, weights)}
    
    GHG_array = np.array(GHG)
    country_GHG_dict = {country_name: GHG for country_name, GHG in zip(country_names, GHG)}
    
    # Final arrays returned to the model: same order as country_weight_dict/
    # country_GHG_dict, i.e. aligned with country_names
    weights_array = np.array(list(country_weight_dict.values()))
    #print("Weights array:", weights_array)
    
    GHG_array = np.array(list(country_GHG_dict.values()))
    # Print the dictionary
    #print("Country-Weight Dictionary:", country_weight_dict)

    return weights_array, GHG_array
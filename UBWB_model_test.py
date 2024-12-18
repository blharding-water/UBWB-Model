# -*- coding: utf-8 -*-
"""
Test function for HD_model.py

Created on Tue Nov 12 17:55:50 2024

@author: bhard
License: CC BY-SA 4.0, https://creativecommons.org/licenses/by-sa/4.0/
"""
import pandas as pd
from UBWB_model import simulate_trace
import UBWB_model_utilities as Mu

def main(data_path,output_path):
    # Codes to validate the 2007HD results and to validate against previous
    # outputs to catch introduced bugs
    
    # Validation against 2007HD, runs 2 and 6, 1930-1963
    # 1930 starts with reservoirs full. Shortages occur in 1963 and 1964.
    # 33 years of accumulation will reveal any differences.
    # 1977 shortages/curtailments will not match due to difference in 
    #   conceptual model. 
    HD_flows = pd.DataFrame(
        columns=["year", "flow"],
        data=[
            [1929, 21829585], [1930, 14621041], [1931, 8474134],
            [1932, 17422187], [1933, 12183500], [1934, 6178192],
            [1935, 12630349], [1936, 14648873], [1937, 14306056],
            [1938, 18148319], [1939, 11164059], [1940, 9931657],
            [1941, 20116678], [1942, 17225136], [1943, 13731401],
            [1944, 15369422], [1945, 14140528], [1946, 11095453],
            [1947, 16439486], [1948, 15139294], [1949, 16933584],
            [1950, 13140416], [1951, 12505894], [1952, 20805422],
            [1953, 11165419], [1954, 8496102], [1955, 9413908],
            [1956, 11426874], [1957, 21500963], [1958, 15862511],
            [1959, 9598169], [1960, 11524160], [1961, 10010259],
            [1962, 17377609], [1963, 8840900], [1964, 10863586],
            [1965, 19875027],[1966,10679844],[1967,11670830],
            [1968,13739932], [1969,15272159], [1970,15344136],
            [1971,15493659], [1972,13186637], [1973,18650193],
            [1974,13285426], [1975,17072661], [1976,11313561],
            [1977,5551188], [1978,15335909], [1979,17825429],
            [1980,17927076], [1981,9015200], [1982,17489400],
            [1983,24361989], [1984,25359376], [1985,21246109],
            [1986,23013446], [1987,15640478], [1988,11456357],
            [1989,9921847], [1990,9639803], [1991,12170021],
            [1992,10895580], [1993,18160118], [1994,11125503],
            [1995,20047166], [1996,14502293], [1997,21622438],
            [1998,16798378], [1999,15934210], [2000,10646526]
        ]
    )
    
    #**********************Run 2********************************

    HD_2007_run2_validation_outputs = simulate_trace(
        HD_flows,
        res_model='active',
        lees_ferry_ann_q=8250000,  # Consistency with 2007HD
        ub_demand=5790000,
        ppr_volume=0
    )
    
    trial_curtailment = HD_2007_run2_validation_outputs[
                        HD_2007_run2_validation_outputs["year"] == 1963
                        ].iloc[0]["curtailment"]
    true_curtailment = 1153349.0 # From 2007HD
    delta = round((trial_curtailment - true_curtailment)
                  / true_curtailment * 100, 3)
    
    print('***HD 2007 Run 2 validation:')
    print(
        f'Curtailment validation:\n'
        f'1963 trial: {trial_curtailment} '
        f' true: {true_curtailment} '
        f'delta: {delta}%'
    )
    trial_curtailment = HD_2007_run2_validation_outputs[
                        HD_2007_run2_validation_outputs["year"] == 1977
                        ].iloc[0]["curtailment"]
    true_curtailment = 3136608  # From 2007HD
    delta = round((trial_curtailment - true_curtailment)
                  / true_curtailment * 100, 3)
    
    print(
        f'1977 trial: {trial_curtailment}'
        f' true: {true_curtailment} '
        f'delta: {delta}%'
    )
    print(
        f'HD 2007 Run 2 mass balance: '
        f'{Mu.check_mass_balance(HD_2007_run2_validation_outputs)}'
    )
    run_name = "2007HD Run 2 validation"
    metadata = ('reservoir_capacity,active,\n'
                'Lees_Ferry_Ann_Q,8250000, MOR\n'
                'UB_demand,5790000\nPPR_volume,0\n'
                'Trigger,False'
                )
    Mu.process_single_trace(
        HD_2007_run2_validation_outputs,
        f'{run_name}',
        output_path,
        metadata = metadata)
    
    #*******************************Run 6********************************
    HD_2007_run6_validation_outputs = simulate_trace(
        HD_flows,
        res_model='live',
        lees_ferry_ann_q=8250000,  # Consistency with 2007HD
        ub_demand=5980000,
        ppr_volume=0
    )
    
    trial_curtailment = HD_2007_run6_validation_outputs[
                        HD_2007_run6_validation_outputs["year"] == 1963
                        ].iloc[0]["curtailment"]
    true_curtailment = 703237 # From 2007HD
    delta = round((trial_curtailment - true_curtailment)
                  / true_curtailment * 100, 3)
    
    print('\n***HD 2007 Run 6 validation:')
    print(
        f'Curtailment validation:\n'
        f'1963 trial: {trial_curtailment} '
        f' true: {true_curtailment} '
        f'delta: {delta}%'
    )
    trial_curtailment = HD_2007_run6_validation_outputs[
                        HD_2007_run6_validation_outputs["year"] == 1977
                        ].iloc[0]["curtailment"]
    true_curtailment = 3665093   # From 2007HD
    delta = round((trial_curtailment - true_curtailment)
                  / true_curtailment * 100, 3)
    
    print(
        f'1977 trial: {trial_curtailment}'
        f' true: {true_curtailment} '
        f'delta: {delta}%'
    )
    print(
        f'HD 2007 Run 6 mass balance: '
        f'{Mu.check_mass_balance(HD_2007_run6_validation_outputs)}'
    )
    run_name = "2007HD Run 6 validation"
    metadata = ('reservoir_capacity,active,\n'
                'Lees_Ferry_Ann_Q,8250000, MOR\n'
                'UB_demand,5980000\nPPR_volume,0\n'
                'Trigger,False'
                )
    Mu.process_single_trace(
        HD_2007_run6_validation_outputs,
        f'{run_name}',
        output_path,
        metadata = metadata)
    #********************Run extreme low flow test******************
    # test inflow less than depletions
    test_flows = pd.DataFrame(columns=["year","flow"],
        #Meko et al., 2007 flows
        data = [
        #year, flow (af)
        [1869,15940000], [1870,12800000], [1871,8560000],
        [1872,16380000], [1873,4000000],  [1874,11660000],
        [1875,13150000], [1876,15120000], [1877,13110000],
        [1878,12710000], [1879,4000000],  [1880,13610000],
        [1881,12330000], [1882,10010000], [1883,11670000],
        [1884,17930000], [1885,17840000], [1886,14150000],
        [1887,9180000], [1888,13940000], [1889,12790000],
        [1890,15430000], [1891,16090000]
        ])
    
    low_flow_test_outputs =   simulate_trace(
                                   test_flows,
                                   res_model = 'active',
                                   lees_ferry_ann_q = 8230000,
                                   ub_demand = 5790000,
                                   ppr_volume = 2267000)
    run_name = "low_flow_test"
    metadata = ('reservoir_capacity,active,\n'
                'Lees_Ferry_Ann_Q,8230000\n'
                'UB_demand,5790000\nPPR_volume,3317000\n'
                'Trigger,False'
                )
    
    Mu.process_single_trace(low_flow_test_outputs,run_name,output_path,metadata = metadata)
    
    print('\n****** Low-flow test *******')
    trial_10yr = low_flow_test_outputs[
                 low_flow_test_outputs["year"] == 1882
                 ].iloc[0]["LF_10yr_flows"]
    
    true_10yr = 81680120  # taken from output of last code revision
    delta = round((trial_10yr - true_10yr) / true_10yr * 100, 3)
    print(
        f"10yr validation: trial: {trial_10yr} true: {true_10yr} "
        f"delta: {delta}%"
    )
    old_mass_balance = 6 # Hand calculated from output of last code revision.
    print(
        f'Low flow test mass balance: trial '
        f'{Mu.check_mass_balance(low_flow_test_outputs)} '
        f'last: {old_mass_balance}'
        )

    # ***********************Test using Meko outputs*******************

    data_file = "meko_et_al_2007_762_2005_trace.csv"
    Meko_LFflows = pd.read_csv(f"{data_path}{data_file}", comment="#")
    Meko_LFflows = Meko_LFflows.astype(int)
    
    Meko_validation_outputs = simulate_trace(
        Meko_LFflows, res_model = 'active',
        lees_ferry_ann_q=8230000,
        ub_demand=5790000,
        ppr_volume=2317000
    )
    
    print('\n***Meko 2007 test:')
    trial_ppr = Meko_validation_outputs[
                Meko_validation_outputs["year"] == 1902
                ].iloc[0]["UB_BU"]
    
    true_ppr = 2317000  # will equal argument passed to simulate_trace
    delta = round((trial_ppr - true_ppr) / true_ppr * 100, 3)
    print(
        f"PPR validation: trial: {trial_ppr} true: {true_ppr} "
        f"delta: {delta}%"
    )
    
    trial_10yr = Meko_validation_outputs[
                 Meko_validation_outputs["year"] == 1902
                 ].iloc[0]["LF_10yr_flows"]
    
    true_10yr = 77847390  # taken from output of last code revision
    delta = round((trial_10yr - true_10yr) / true_10yr * 100, 3)
    print(
          f"10yr validation: trial: {trial_10yr} true: {true_10yr} "
          f"delta: {delta}%"
          )
    
    old_mass_balance = -57 # Hand calculated from output of last code revision.
    print(
        f'Meko 2007 mass balance: trial '
        f'{Mu.check_mass_balance(Meko_validation_outputs)}'
        f' last: {old_mass_balance}'
    )
    run_name = "Meko 2007 paleo_5790_8230_2317_PP_MOR"
    metadata = ('reservoir_capacity,active,\n'
                'Lees_Ferry_Ann_Q,8230000, MOR\n'
                'UB_demand,5790000\nPPR_volume,3317000\n'
                'Trigger,False'
                )
    Mu.process_single_trace(
        Meko_validation_outputs,
        f'{run_name}',
        output_path,
        metadata = metadata) 

if __name__ == '__main__':

    data_path = './'
    output_path = './'
    
    main(data_path,output_path)

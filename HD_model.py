import pandas as pd
import numpy as np
import analyze_spells

TOLERANCE = 5 # acre-feet criterion for closure of evaporation solution
MAX_TRIALS = 5


def evaporation(reservoir_contents):
    """
    Return reservoir evaporation as a function of reservoir contents.
    Developed using regression against HD 2007 outputs. This equation 
    differs from equations called out in HD 2007 text.  USBR seems to base
    evap on total storage, even though they say they base it on CRSP 
    storage. This equation results from regression against total storage,
    using the average of start & end reservoir contents.
    """
    # This is the equation to use if simulating use of power pools:
    #    return(int(round(0.021292*reservoir_contents+5016.897167,0)))
    # If power pools are to be protected use this equation: 
    return int(round(0.020874 * reservoir_contents + 132877, 0))


def mor_release(lf_ann_q, lf_deficit, *args):
    """
    Minimum Objective Release as in LROC
    but values other than 8.23 can be provided.
    """ 
    return max(lf_ann_q, lf_deficit)


def no_mor_release(lf_ann_q, lf_deficit, *args):
    """
    Release the compact deficit.
    """
    return lf_deficit


def trigger_cutback(res_capacity, res_contents, ub_non_ppr_depls):
    """
    An arbitrary trigger scheme based on reservoir contents
    used for sensitivity analysis of the efficacy of triggers.
    
    Returns the amount to cut back Upper Basin beneficial use.
    """
    state = res_contents / float(res_capacity)
    if state > 0.33:
        return 0
    elif state > 0.25:
        return 0.1 * ub_non_ppr_depls
    elif state > 0.15:
        return 0.2 * ub_non_ppr_depls
    elif state > 0.1:
        return 0.4 * ub_non_ppr_depls
    else:
        return 0.7 * ub_non_ppr_depls


def check_mass_balance(run_output):
    """
    Calculate mass balance around a run

    Parameters
    ----------
    run_output : pandas.Dataframe
        Output from simulate_trace.

    Returns
    -------
    mass_balance : same as values in output
        DESCRIPTION.
    """
    
    mass_balance = run_output.iloc[0].start_con - run_output.iloc[-1].end_con
    mass_balance += run_output['inflow'].sum()
    mass_balance -= run_output['UB_BU'].sum()
    mass_balance -= run_output['evap'].sum()
    mass_balance -= run_output['LF_flow'].sum()
    return mass_balance


def simulate_trace(input_data, start_contents, reservoir_capacity=29530030,
                   lees_ferry_ann_q=8230000, nyrs=10,
                   lees_ferry_n_year_record=None,
                   lf_release=mor_release, ub_demand=5790000,
                   ppr_volume=2267000, trigger_func=None):
    """
    Simulate water flow in the Upper Basin.
    """
    if lees_ferry_n_year_record is None:
        lees_ferry_n_year_record = nyrs * [lees_ferry_ann_q]
    lees_ferry_cum_q = nyrs * lees_ferry_ann_q

    if 'UB demand' not in input_data.columns:
        input_data['UB demand'] = ub_demand

    outputs = pd.DataFrame(columns=[
        "year", "inflow", "start_con", 'trgr_cut', "UB_dmd", "evap",
        "net_avail", "spill", "curtailment", "end_con", "UB_BU", "UB_CU",
        f"LF_{nyrs}yr_flows", "LF_deficit", "LF_flow", "time_from_reset",
        "evap_trials"
    ])

    time_from_reset = 0

    for index, row in input_data.iterrows():
        inflow = row["flow"]
        year = row["year"]
        outputs.loc[index, "year"] = year
        outputs.loc[index, "inflow"] = inflow
        outputs.loc[index, "start_con"] = start_contents

        ub_depletions = row['UB demand']
        if trigger_func:
            cutback = trigger_func(reservoir_capacity, start_contents,
                                   ub_depletions - ppr_volume)
            outputs.loc[index, 'trgr_cut'] = cutback
            ub_depletions -= cutback
        # For the case where inflows are less than Upper Basin demand
        ub_depletions = min(inflow, ub_depletions)
        outputs.loc[index, "UB_dmd"] = ub_depletions
        
        # Look back nyrs-1 years
        lees_ferry_n_year_record.pop()  
        # in order to calculate this year's flow requirement
        lees_ferry_deficit = max(0, 
                                (lees_ferry_cum_q 
                                 - sum(lees_ferry_n_year_record))
                                )
        
        outputs.loc[index, "LF_deficit"] = lees_ferry_deficit
        lf_target = lf_release(lees_ferry_ann_q, lees_ferry_deficit)

        depleted_inflow = inflow - ub_depletions
        evap = evaporation(start_contents)
        evap_trial = 0

        while True:
            evap_trial += 1
            trial_evap = evap
            available_to_store = (depleted_inflow
                                  + start_contents
                                  - lf_target
                                  - trial_evap)
            
            trial_contents = max(
                                 min(available_to_store, reservoir_capacity),
                                 0)
            
            evap = evaporation((start_contents + trial_contents) / 2)

            if abs(trial_evap - evap) < TOLERANCE or evap_trial > MAX_TRIALS:
                break

        end_contents = trial_contents
        outputs.loc[index, "evap_trials"] = evap_trial
        outputs.loc[index, "evap"] = evap
        outputs.loc[index, "end_con"] = end_contents
        outputs.loc[index, "net_avail"] = available_to_store

        spill = max(available_to_store - reservoir_capacity, 0)
        outputs.loc[index, "spill"] = spill
        
        # Address case where inflows are less than PPR volume
        ppr_supply = min(ppr_volume, inflow - evap)
        
        # Curtail by the shortfall meeting LF target
        # or the amount of non-PPR use, whichever is less
        curtailment = min(-min(available_to_store, 0),
                          ub_depletions 
                          - min(ub_depletions, ppr_supply))
        
        outputs.loc[index, "curtailment"] = curtailment
        
        # Calculate the water balance to determine the flow at Lee Ferry,
        # which is the release from the reservoir
        lees_ferry_flow = int(round(depleted_inflow
                                    + start_contents 
                                    - end_contents 
                                    - evap 
                                    + curtailment, -1))
        
        lees_ferry_n_year_record.insert(0, lees_ferry_flow)
        outputs.loc[index, "LF_flow"] = lees_ferry_flow
        outputs.loc[index, f"LF_{nyrs}yr_flows"] = sum(lees_ferry_n_year_record)

        ub_bu = int(round(ub_depletions - curtailment, 0))
        ub_cu = int(round(ub_bu + evap, 0))
        outputs.loc[index, "UB_BU"] = ub_bu
        outputs.loc[index, "UB_CU"] = ub_cu
           
        # If we have a curtailment, how long has it been from the last 
        # spill or curtailment?
        if (int(round(curtailment))!= 0):
            if time_from_reset > 0:
                outputs.at[index,"time_from_reset"] = time_from_reset
            time_from_reset = 0
        elif (spill != 0):
            time_from_reset = 0
        else:
            time_from_reset += 1            

        start_contents = end_contents

    return outputs
    

def calculate_intervals(data):
    """
    Calculate intervals based on spill and curtailment events.
    
    The interval includes the year of the curtailment since depletions
    accumulate in that year.
    """
    events = pd.DataFrame(
        columns=["year", "spill", "curtailment", 
                 "from_spill", "from_curtailment"]
    )
    
    from_spill = np.nan
    from_curtailment = np.nan
    s_intervals = []
    c_intervals = []
    
    row_no = 0
    for _, row in data.iterrows():
        spill = row["spill"]
        curtailment = row["curtailment"]

        if spill > 0 or curtailment > 0:
            events.at[row_no, "year"] = row["year"]
            events.at[row_no, "spill"] = spill
            events.at[row_no, "curtailment"] = curtailment

            if curtailment > 0:
                events.at[row_no, "from_spill"] = from_spill
                events.at[row_no, "from_curtailment"] = from_curtailment
                s_intervals.append(from_spill)
                c_intervals.append(from_curtailment)
                from_curtailment = 1
                from_spill = np.nan
            
            if spill > 0:
                from_spill = 1
                from_curtailment = np.nan
        else:
            from_spill += 1 if not pd.isna(from_spill) else 0
            from_curtailment += 1 if not pd.isna(from_curtailment) else 0
        
        row_no += 1

    return events, s_intervals, c_intervals

def process_single_trace(outputs, run_name, output_path, metadata=''):
    """
    Process a single trace of outputs, saving time series and curtailment data.
    """
    # Prepare metadata
    metadata = f"{run_name}\n{metadata}\n"
    
    # Write time series CSV
    ts_file_path = f"{output_path}{run_name}.TS.csv"
    with open(ts_file_path, 'w', newline='\n') as f:
        f.write(metadata)
        outputs.to_csv(f, index=False, lineterminator='')
    
    # Calculate intervals
    events, s_intervals, c_intervals = calculate_intervals(outputs)

    # Analyze spells
    nested_spells = {}
    spells = {}
    analyze_spells.characterize_spells(
        outputs['curtailment'], nested_spells, spells
    )
    
    # Write curtailment data
    curtailment_file_path = f"{output_path}{run_name}.curtailments.csv"
    with open(curtailment_file_path, 'w', newline='\n') as outfile:
        outfile.write(f"Outputs for {run_name}\n")
        
        # Write spell information
        outfile.write(
            f"Spell information for {run_name}\n"
            "Discrete spell events\n"
        )
        analyze_spells.write_spell_dict(outfile, spells)
        outfile.write("Nested spells\n")
        analyze_spells.write_spell_dict(outfile, nested_spells)
        
        # Write curtailment intervals
        outfile.write(
            "Curtailment intervals from spill/curtailment\n"
        )
        
        if not s_intervals:
            outfile.write("Spill: Max: ,--, Min: ,--, Mean: ,--\n")
        else:
            outfile.write(
                f"Spill: Max:,{np.nanmax(s_intervals)}, "
                f"Min:,{np.nanmin(s_intervals)}, "
                f"Mean:,{np.nanmean(s_intervals):.2f}\n"
            )
        
        if not c_intervals:
            outfile.write("Curt.: Max: ,--, Min: ,--, Mean: ,--\n")
        else:
            outfile.write(
                f"Curt.: Max:,{np.nanmax(c_intervals)}, "
                f"Min:,{np.nanmin(c_intervals)}, "
                f"Mean:,{np.nanmean(c_intervals):.2f}\n"
            )



def main(data_path,output_path):
    # Codes to validate the 2007HD results and to validate against prefious
    # outputs to catch introduced bugs
    
    # Validation against 2007HD, Run 2, 1930-1963
    # 1930 starts with reservoirs full. Shortages occur in 1963 and 1964.
    # 33 years of accumulation will reveal any differences.
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
            [1965, 19875027]
        ]
    )

    HD_2007_validation_outputs = simulate_trace(
        HD_flows, 29530030,
        reservoir_capacity=29530030,
        lees_ferry_ann_q=8250000,  # Consistency with 2007HD
        ub_demand=5790000,
        ppr_volume=0
    )
    
    trial_curtailment = HD_2007_validation_outputs[
                        HD_2007_validation_outputs["year"] == 1963
                        ].iloc[0]["curtailment"]
    true_curtailment = 1153349.0 # From 2007HD
    delta = round((trial_curtailment - true_curtailment)
                  / true_curtailment * 100, 3)
    
    print('***HD 2007 validation:')
    print(
        f'Curtailment validation: trial: {trial_curtailment} '
        f'true: {true_curtailment} '
        f'delta: {delta}%'
    )
    print(
        f'HD 2007 mass balance: '
        f'{check_mass_balance(HD_2007_validation_outputs)}'
    )

    # Validate using Meko outputs

    data_file = "meko_et_al_2007_762_2005_trace.csv"
    Meko_LFflows = pd.read_csv(f"{data_path}{data_file}", comment="#")
    Meko_LFflows = Meko_LFflows.astype(int)
    
    Meko_validation_outputs = simulate_trace(
        Meko_LFflows, 29530030,
        reservoir_capacity=29530030,
        lees_ferry_ann_q=8230000,
        ub_demand=5790000,
        ppr_volume=2317000
    )
    
    print('***Meko 2007 validation:')
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
        f'{check_mass_balance(Meko_validation_outputs)}'
        f' last: {old_mass_balance}'
    )
    Meko_validation_outputs.to_csv(
        f"{output_path}meko_2007_validation_outputs.csv"
        )
    run_name = "Meko 2007 paleo_5790_8230_2317_PP_MOR"
    metadata = ('reservoir_capacity,29530030,\n'
                'Lees_Ferry_Ann_Q,8230000\n'
                'UB_demand,5790000\nPPR_volume,3317000\n'
                'Trigger,False'
                )
    process_single_trace(
        Meko_validation_outputs,
        f'{run_name}.test',
        output_path,
        metadata = metadata)
    
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
                                   29530030,
                                   reservoir_capacity = 29530030,
                                   lees_ferry_ann_q = 8230000,
                                   ub_demand = 5790000,
                                   ppr_volume = 2267000)
    run_name = "low_flow_test"
    metadata = ('reservoir_capacity,29530030,\n'
                'Lees_Ferry_Ann_Q,8230000\n'
                'UB_demand,5790000\nPPR_volume,3317000\n'
                'Trigger,False'
                )
    
    process_single_trace(low_flow_test_outputs,run_name,output_path,metadata = metadata)
    
    print('****** Low-flow test *******')
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
        f'{check_mass_balance(low_flow_test_outputs)} '
        f'last: {old_mass_balance}'
        )
    
if __name__ == '__main__':
    
    output_path = (
        "C:/Users/bhard/BLH Folders/technical/computer codes/"
        "_my_python_library/HD_model working/"
    )
    
    data_path = (
        "C:/Users/bhard/BLH Folders/technical/Colorado River/"
        "2019 Risk Analysis/HD model 2019/inflows/lees_ferry_(gage)/"
        )
    
    main(data_path,output_path)


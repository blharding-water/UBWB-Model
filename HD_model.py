"""
Python 3
NOTE: This code has errors at boundary cases and certainly has undetected errors
	  Use at your own risk.
	  
Title: Colorado River HD Model
Author: Ben Harding, bharding@lynker.com
Source:
License: CC BY-SA 4.0, https://creativecommons.org/licenses/by-sa/4.0/

Program to simulate the Upper Basin of the Colorado River using the approach
of the 2007 USBR Hydrologic Determination method (HD2007).

 Author: Ben Harding, 2008
 1-14-2010 CRWAS ensemble method
 9-22-2019 Single-trace method
 9-26-2019 refactoring:    1 make trace method module w/ integral validataion
                           2 use pandas
 4-25-2024 Added variable for accumulation period for Lee Ferry non-depletion
           obligation.  
 6-17-2024 Added code to simulate a trigger system
 6-20-2024 Added code to use 'UB demand' time series if present in input_data
 7-2-2024 Added code to eliminate negative shortages when demand is set
          less than PPR amount. This has happened in an automated yield 
          determination but it should not happen in a real case.But, JIC.
 8-11-2024 Added code to check mass balance
 10-29-2024 Added code to reduce UB_depletions when inflows are insufficient
            This could happen in a climate-change case.
 11-12-2024 Refactored into three modules, HD_model.py, HD_model_utilties.py
            and HD_model_test.py

"""

import pandas as pd
import HD_model_utilities as HDu

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
    

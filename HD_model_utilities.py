# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 17:25:44 2024

@author: Ben Harding

"""
import pandas as pd
import numpy as np
import sys

"""
Output Utilities
"""

# Define default quantiles
DEFAULT_QUANTILES = [10, 25, 50, 75, 90]

def write_percentiles(file, data, labels, quantiles=DEFAULT_QUANTILES):
    """Write percentiles, max, min, and mean values of data to a CSV file."""
    # Data is a tuple of lists
    columns = quantiles + ['Max', 'Min', 'Mean']
    out_df = pd.DataFrame(columns=columns)
    
    for i, series in enumerate(data):
        # Calculate percentiles and summary statistics
        percentiles = list(np.nanpercentile(series, quantiles))
        percentiles.append(np.nanmax(series))
        percentiles.append(np.nanmin(series))
        percentiles.append(float(np.nanmean(series)))
        
        # Assign calculated values to the DataFrame row
        out_df.loc[labels[i]] = percentiles
    
    out_df.to_csv(file, lineterminator="\n")
    return out_df


def write_ts_output(file_spec, run_metadata, outputs):
    """Write metadata and time series outputs to a CSV file."""
    with open(file_spec, 'w') as file:
        for key, value in run_metadata.items():
            file.write(f"{key},{value}\n")

    with open(file_spec, 'a') as file:
        outputs.to_csv(file, index=False, lineterminator='\n')

'''
Spell utilities
'''
def characterize_spells(data, n_spell_dict=None, i_spell_dict=None):
    """
    Function to analyze spells of events in a series.
    
    Arguments:
        data: list or Pandas series containing an ordered sequence
              of positive excursions and zeros.
        n_spell_dict (opt): dictionary holding lists of nested spells keyed by duration.
        i_spell_dict (opt): dictionary holding lists of independent spells keyed by duration.
    
    Returns:
        n_spell_dict, i_spell_dict

    Events are positive excursions from zero, so for conventional
    analysis of precipitation or streamflow spells the time series
    must be normalized by subtracting the time-series values from the 
    time-series mean.

    Spells are reported by the average value of the excursions over the 
    duration of the spell.

    Two types of spells are characterized: conventional independent spells
    and 'nested' spells (spells of all durations contained within independent spells).
    """
    if n_spell_dict is None:
        n_spell_dict = {}
    if i_spell_dict is None:
        i_spell_dict = {}

    remainder = list(data)  # will accept a series
    while True:
        # Strip off leading zeros
        while remainder and remainder[0] == 0:
            remainder.pop(0)
        if not remainder:
            break

        # Extract and catalog the independent spell
        try:
            next_zero = remainder.index(0)
        except ValueError:
            next_zero = len(remainder)

        run = remainder[:next_zero]
        duration = len(run)
        if duration not in i_spell_dict:
            i_spell_dict[duration] = []
        i_spell_dict[duration].append(sum(run) / float(duration))

        # Find and catalog the nested spells
        for nested_duration in range(1, len(run) + 1):
            for idx in range(0, len(run) - nested_duration + 1):
                if nested_duration not in n_spell_dict:
                    n_spell_dict[nested_duration] = []
                n_spell_dict[nested_duration].append(
                    sum(run[idx:idx + nested_duration]) / float(nested_duration)
                )

        remainder = remainder[next_zero:]

    return n_spell_dict, i_spell_dict


def write_spell_dict(file, spell_dict, title=None):
    """Writes each duration on one line."""
    if title:
        file.write(f"{title}\n")
    file.write('Duration\n')
    for key, values in spell_dict.items():
        values_str = ",".join(map(str, values))
        file.write(f"{key}, {values_str}\n")


def write_spell_percentiles(file, spell_dict, quantiles, title=None):
    """Writes percentiles for each spell duration."""
    labels = list(spell_dict.keys())
    data = list(spell_dict.values())
    if title:
        file.write(f"{title}\n")
    write_percentiles(file, data, labels, quantiles)

'''
HD model functions
'''

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
    characterize_spells(
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
        write_spell_dict(outfile, spells)
        outfile.write("Nested spells\n")
        write_spell_dict(outfile, nested_spells)
        
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
            
'''
Test functions for this module\.
'''       
        
def test_output_utilities(output_path):
    # Sample data and labels for testing
    data = (range(0, 1001), range(1, 10))
    labels = (1, 2)
    
    # Test writing percentiles to standard output
    write_percentiles(sys.stdout, data, labels)
    
    # Specify output path (change as needed)
    write_percentiles(f"{output_path}test_percentiles.csv", data, labels)
    
    # Additional test with a subset of data
    d = (data[1],)
    write_percentiles(sys.stdout, d, labels)       
    
def test_spell_utilities(output_path):
    
    test_data = [0, 0, 1, 3, 0, 0, 4, 5, 6, 0, 0, 7, 8, 9, 10,
                0, 0, 11, 12, 13, 14, 15, 0]

    nested_spell_dict, independent_spell_dict = characterize_spells(test_data)

    print(f'len test data: {len(test_data)}.')
    print(f'data: {test_data}')
    print('test data results:')
    print('****nested spell dict')
    print(f'correct: 2: {[2.0, 4.5, 5.5, 7.5, 8.5, 9.5, 11.5, 12.5, 12.5, 13.5, 14.5]}')
    write_spell_dict(sys.stdout, nested_spell_dict, title='nested spells')

    print('\n****independent spell dict')
    print('correct: {2: [2.0], 3: [5.0], 4: [8.5], 5: [13.0]}')
    write_spell_dict(sys.stdout, independent_spell_dict, title='independent spells')

    with open(output_path + "test_write_spell_percentiles.csv", "w") as file:
        quantiles = [10, 25, 50, 75, 90]
        write_spell_percentiles(file, nested_spell_dict, quantiles, title='nested spell percentiles')
    
if __name__ == '__main__':

    output_path = './'
    test_output_utilities(output_path)
    test_spell_utilities(output_path)
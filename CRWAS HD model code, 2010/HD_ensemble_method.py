# program to simulate USBR Hydrologic Determination method against ensemble of flow traces
# For CRWAS
import sys, text_data_files
from build_traces import file_names

TRUE = 1==1
FALSE = not TRUE

hydrology_directory_path = r"P:/589 CRWAS/7.12/AnnualFlowWY/"
HD_path = "P:/589 CRWAS/8.3 Big River Model/HD Method/"
historical_hydrology_file_name = "LeesFerryNFWY1950-2005.txt"
historical_hydrology_file_name = "LeeFerryNFWY1950-2005.txt"
full_historical_hydrology_file_name = "HD_natural_flows_list.txt"
full_historical_hydrology_file_name = "HD_Lee_Ferry_natural_flows_list.txt"
extended_historical_hydrology_file_name = "LeesFerryNFWY1950-2005_ensemble.txt"
extended_as_if_hydrology_file_names = {2040:("sresa1b.ncar_ccsm3_0.2_2040.Adj_allStat.AcFt.txt_WY_AnnualTotals.SumToLeesFerry_ensemble.txt",
                                    "sresa2.miroc3_2_medres.1_2040.Adj_allStat.AcFt.txt_WY_AnnualTotals.SumToLeesFerry_ensemble.txt",
                                    "sresa2.mri_cgcm2_3_2a.1_2040.Adj_allStat.AcFt.txt_WY_AnnualTotals.SumToLeesFerry_ensemble.txt",
                                    "sresa2.ncar_pcm1.3_2040.Adj_allStat.AcFt.txt_WY_AnnualTotals.SumToLeesFerry_ensemble.txt",
                                    "sresb1.cccma_cgcm3_1.2_2040.Adj_allStat.AcFt.txt_WY_AnnualTotals.SumToLeesFerry_ensemble.txt"),
                              2070:("sresa1b.gfdl_cm2_0.1_2070.Adj_allStat.AcFt.txt_WY_AnnualTotals.SumToLeesFerry_ensemble.txt",
                                    "sresa1b.mri_cgcm2_3_2a.4_2070.Adj_allStat.AcFt.txt_WY_AnnualTotals.SumToLeesFerry_ensemble.txt",
                                    "sresa1b.ncar_ccsm3_0.2_2070.Adj_allStat.AcFt.txt_WY_AnnualTotals.SumToLeesFerry_ensemble.txt",
                                    "sresa2.ncar_pcm1.3_2070.Adj_allStat.AcFt.txt_WY_AnnualTotals.SumToLeesFerry_ensemble.txt",
                                    "sresb1.mpi_echam5.1_2070.Adj_allStat.AcFt.txt_WY_AnnualTotals.SumToLeesFerry_ensemble.txt")
                              }
tolerance = .1
max_trials = 5

# Simulation system parameters
Lees_Ferry_Annual_Delivery = 8250000
##Lees_Ferry_Annual_Delivery = 7500000
UB_depletions = 5980000  # HD at 82.5
##UB_depletions = 6760000  # HD at 75.0
##UB_depletions = 5379000    # max 2060 requests
PPR_volume = 0   # may change as a sensitivity analysis
reservoir_capacity = 33833590  # all reservoirs
##reservoir_capacity = 24322000   # Powell only
CRSP_fraction = 0.861607148
##CRSP_fraction = 1.0             # use with powell-only capacity.  All evap happens there


# initialize system variables
Lees_Ferry_cumulative_flow_requirement = 10* Lees_Ferry_Annual_Delivery
Lees_Ferry_10_year_record = 10*[Lees_Ferry_Annual_Delivery,]

def evaporation(reservoir_contents):
    return 0.0247*reservoir_contents*CRSP_fraction+5017

def stopping_rule(start_values,end_values,trial):
##    print "testing stopping rule:", start_values,end_values,trial
    if trial > max_trials:
        return TRUE
    if (abs(start_values[0]-end_values[0])/reservoir_capacity <= tolerance and
       abs(start_values[1]-end_values[1])/Lees_Ferry_cumulative_flow_requirement <= tolerance):
        return TRUE


def simulate_trace(flow_data):
    global reservoir_contents, inflow
    reservoir_contents = reservoir_capacity
    trial = 0
    start_values = (reservoir_contents,sum(Lees_Ferry_10_year_record))
    while TRUE:
        trial += 1
        curtailments = []
        cumulative_inflow = cumulative_outflow = 0
        for inflow in flow_data:
            cumulative_inflow += inflow
            available_to_store = (inflow+reservoir_contents
                                  -UB_depletions-Lees_Ferry_Annual_Delivery
                                  -evaporation(reservoir_contents)
                                  )
            spill = max(available_to_store - reservoir_capacity, 0.0)
            shortage = min(available_to_store, 0.0)
            Lees_Ferry_flow = Lees_Ferry_Annual_Delivery+spill+shortage # only one will be non-zero
            Lees_Ferry_10_year_record.pop()   # get rid of 11th year
            trial_sum = sum(Lees_Ferry_10_year_record) + Lees_Ferry_flow
            if trial_sum < Lees_Ferry_cumulative_flow_requirement:
                deficit = Lees_Ferry_cumulative_flow_requirement - trial_sum
                curtailment = min(deficit,UB_depletions-PPR_volume,inflow-evaporation(reservoir_contents))
            else:
                curtailment = 0
            Lees_Ferry_flow = Lees_Ferry_flow + curtailment
            cumulative_outflow += (UB_depletions-curtailment+Lees_Ferry_flow+evaporation(reservoir_contents))
            curtailments.append(curtailment)
            Lees_Ferry_10_year_record.insert(0,Lees_Ferry_flow)
            Lees_Ferry_10_year_cumulative_flow = sum(Lees_Ferry_10_year_record)
##            print trial,inflow,reservoir_contents,evaporation(reservoir_contents),available_to_store,spill,shortage,trial_sum,curtailment,Lees_Ferry_flow,Lees_Ferry_10_year_cumulative_flow,Lees_Ferry_10_year_record
            reservoir_contents = max(min(available_to_store, reservoir_capacity),0)  # re-set for next run    
        end_values = (reservoir_contents,sum(Lees_Ferry_10_year_record))
        if stopping_rule(start_values,end_values,trial):
            cumulative_inflow += (start_values[0] - end_values[0])  # reservoir change in storage
            return curtailments,sum(curtailments)/len(curtailments),trial,cumulative_inflow,cumulative_outflow,
        start_values = end_values
            
    
# run full historical data
flow_data = text_data_files.textdatafile("%s%s" % (HD_path,full_historical_hydrology_file_name)).readdatarecords()
curtailments, avg_curtailment,trials,cumulative_inflow,cumulative_outflow = simulate_trace(flow_data[0])
print "Full historical: average curtailment %s trials %s in %s out %s  " %  (avg_curtailment, trials,cumulative_inflow,cumulative_outflow)
print "\n".join(map(str,curtailments))

# run historical study period, 1950-2005
fd = text_data_files.textdatafile("%s%s" % (hydrology_directory_path,file_names[0][0])).readdatarecords()
flow_data = []
for year in fd:
    flow_data.append(year[1])
curtailments, avg_curtailment,trials,cumulative_inflow,cumulative_outflow = simulate_trace(flow_data)
print "1950-2005 historical: average curtailment %s trials %s in %s out %s " %  (avg_curtailment, trials,cumulative_inflow,cumulative_outflow)
print "\n".join(map(str,curtailments))

### run extendec historical hydrology (100 traces)
##flow_data = text_data_files.textdatafile("%s%s" % (HD_path,extended_historical_hydrology_file_name)).readdatarecords()
##trace_no = 0
##curtailment_population = []
##print "extended historical"
##for trace in flow_data:
##    curtailments, avg_curtailment,trials,cumulative_inflow,cumulative_outflow = simulate_trace(trace)
##    trace_no += 1
##    print trace_no,avg_curtailment,trials,cumulative_inflow,cumulative_outflow
##    curtailment_population.extend(curtailments)

# run all hydrology cases
results_dict = {}
for hydrology_file in file_names:
    file_name = hydrology_file[1]
    run_name = file_name.split("_ensemble.txt")[0]
    flow_data = text_data_files.textdatafile("%s%s" % (HD_path,file_name)).readdatarecords()
    trace_no = 0
    curtailment_population = []
    print "Running %s" % run_name
    for trace in flow_data:
        curtailments, avg_curtailment,trials,cumulative_inflow,cumulative_outflow = simulate_trace(trace)
        trace_no += 1
        print trace_no,avg_curtailment,trials,"%.2f" % (cumulative_inflow-cumulative_outflow)
        curtailment_population.extend(curtailments)
    results_dict[run_name] = curtailment_population

# print curtailments in matrix
outfile = open("%s%s" % (HD_path,"all_curtailment_populations.txt"),"w")
for hydrology_file in file_names:
    outfile.write("%s\t" % hydrology_file[1].split("_ensemble.txt")[0])
outfile.write("\n")
for row in range(len(curtailment_population)):    # use left over to get length
    for hydrology_file in file_names:
        outfile.write("%.3f\t" % results_dict[hydrology_file[1].split("_ensemble.txt")[0]][row])
    outfile.write("\n")
outfile.close()
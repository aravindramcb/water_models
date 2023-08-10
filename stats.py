# Statistical tests to be performed for different results
import os
from collections import defaultdict

import matplotlib.pyplot as plt
from scikit_posthocs import posthoc_dunn
from statsmodels.sandbox.stats.multicomp import MultiComparison

from libs import transport_events_analysis as tt_events
from scipy.stats import kruskal
import pandas as pd
import seaborn as sns


def transit_time_stats(tt_results:str, sim_results:str, gd:dict):
    """
    Process data for histogram and plot transit time per event.
    1 frame = 10ps in my system
    transit time per event = number of frames per event * 10 for a single entry or release event
    :param tt_results: Location of Transport tools results
    :param sim_results: Location of simulation results
    :return: saves a plot of histogram in the specified location
    """
    frames = tt_events.get_transit_time(tt_results=tt_results, simulation_results=sim_results,
                                        groups_definitions=gd,
                                        frame_numbers=False)
    combined = frames[0]
    models = ['opc', 'tip3p', 'tip4pew']
    names = ["1.4", "1.8", "2.4", "3"]
    group_number = [1,2,3,4]
    row, col = (4, 3)  # 4 groups and 3 models

    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.set_context(context="paper", font_scale=1)
    fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(10, 10), dpi=150)

    data_per_group = defaultdict(list)
    data_per_group_split_by_models = []
    for r in range(row):
        data_within_group=defaultdict(list)
        for c in range(col):
            data = []
            current_model = models[c]
            group_name = f"TCG{group_number[r]}"
            # combine data from all sims
            for sim in range(1,6):
                sim_name = f"{names[r]}A_{current_model}_{sim}"
                _sim_time = combined[sim_name]
                data.extend(_sim_time)
            data_within_group[current_model] = data

            # combine opc,tip3p,tip4p data per group for Kruskal-Wallis test
            data_per_group[group_name].extend(data)

        # split data for comparing significance within groups
        data_per_group_split_by_models.append(data_within_group)

    # check significant differences between groups- Kruskal-Wallis test

    # statistical significance betwen groups regardless of models
    statistic, p_value = kruskal(data_per_group["TCG1"],data_per_group["TCG2"],
                                 data_per_group["TCG3"],data_per_group["TCG4"])
    print("P value of Ksuskal-Wallis test regardless of models between groups")
    print(p_value)

    # Statistical significance between groups for different models between groups
    opc_values =[]
    tip3p_values =[]
    tip4pew_values=[]
    for group_dict in data_per_group_split_by_models:
        opc_values.append(group_dict["opc"])
        tip3p_values.append(group_dict["tip3p"])
        tip4pew_values.append(group_dict["tip4pew"])
    statistic, p_value = kruskal(*opc_values)
    print(f"for opc \n p_value={p_value}, statistic={statistic}")
    statistic, p_value = kruskal(*tip3p_values)
    print(f"for tip3p \n p_value={p_value}, statistic={statistic}")
    statistic, p_value = kruskal(*tip4pew_values)
    print(f"for tip4pew \n p_value={p_value}, statistic={statistic}")

    # Posthoc test to compare differences between models within groups

    # Create a list of dictionaries, each representing a key-value pair in the defaultdict
    def perform_test(dict):
        data_list = [{'Group': key, 'Value': value} for key, values in dict.items() for value
                     in values]
        df = pd.DataFrame(data_list)
        result = posthoc_dunn(df, group_col='Group', val_col='Value', p_adjust='bonferroni')
        return result

    group_number = 1
    for group in data_per_group_split_by_models:
        result= perform_test(group)
        group_name =f"TCG{group_number}"
        print(f"For group {group_name} \n -------------------------------------- \n")
        print(result)
        group_number +=1

def percent_events_occurance_statistics(tt_result, sim_results,tunnels_def):
    events = tt_events.fraction_events_occurrence(tt_result, sim_results, tunnels_def)
    combined = events[0:5]
    opc_values = []
    tip3p_values =[]
    tip4pew_values = []
    for sim in combined:
        sim.columns = ["opc","tip3p", "tip4pew"]
        opc_values.append(sim['opc'])
        tip3p_values.append(sim['tip3p'])
        tip4pew_values.append(sim['tip4pew'])

    # to establish statistical difference between groups within models
    statistic, p_value = kruskal(*opc_values)
    print(f"for opc \n p_value={p_value}, statistic={statistic}")
    statistic, p_value = kruskal(*tip3p_values)
    print(f"for tip3p \n p_value={p_value}, statistic={statistic}")
    statistic, p_value = kruskal(*tip4pew_values)
    print(f"for tip4pew \n p_value={p_value}, statistic={statistic}")

    # post hoc test to test within groups
    def perform_test(dict):
        data_list = [{'Group': key, 'Value': value} for key, values in dict.items() for value
                     in values]
        df = pd.DataFrame(data_list)
        result = posthoc_dunn(df, group_col='Group', val_col='Value', p_adjust='bonferroni')
        return result

    group_number = 1
    for group in combined:
        result= perform_test(group)
        group_name =f"TCG{group_number}"
        print(f"For group {group_name} \n -------------------------------------- \n")
        print(result)
        group_number +=1

def average_bottleneck_statistics(csv_file:str):
    bottlenecks= pd.read_csv(csv_file,usecols=[2,3,4,5])

    statistic, p_value = kruskal(bottlenecks["1.4A"], bottlenecks["1.8A"], bottlenecks["2.4A"], bottlenecks["3A"],
                                 nan_policy='omit')

    def perform_test(dict):
        data_list = [{'Group': key, 'Value': value} for key, values in dict.items() for value
                     in values]
        df = pd.DataFrame(data_list)
        result = posthoc_dunn(df, group_col='Group', val_col='Value', p_adjust='bonferroni')
        return result
    result = perform_test(bottlenecks)
    print(result)

def full_bottleneck_statistics():
    # Get the original caver IDs for the given group name (P1,P2,P3) for all simulations
    bottleneck_full = pd.read_csv("/home/aravind/PhD_local/dean/figures/bottlenecks/"
                                  "time_evolution/bottleneck_per_sim.csv")
    save_loc ="/home/aravind/PhD_local/dean/statistics/"
    tcg0=bottleneck_full.iloc[:, 30:45]
    tcg1=bottleneck_full.iloc[:, 0:15]
    tcg2=bottleneck_full.iloc[:, 15:30]
    tcg3= bottleneck_full.iloc[:,45:60]
    tcg4=bottleneck_full.iloc[:,60:75]

    import statsmodels.stats.multicomp as mc
    def perform_test(dataframe):
        col_names = ['opc_1', 'opc_2', 'opc_3', 'opc_4', 'opc_5', 'tip3p_1', 'tip3p_2', 'tip3p_3', 'tip3p_4', 'tip3p_5',
                     'tip4pew_1', 'tip4pew_2', 'tip4pew_3', 'tip4pew_4', 'tip4pew_5']
        dataframe.columns = col_names
        df_long = pd.melt(dataframe.reset_index(drop=True), var_name='Group', value_name='Value')
        dunns_result = posthoc_dunn(df_long,group_col='Group',val_col='Value', p_adjust='bonferroni')
        return dunns_result
    i=0
    for group in [tcg0,tcg1,tcg2,tcg3,tcg4]:
        group_name = f"TCG{i}"
        result= perform_test(group)
        plt.figure(figsize=(10, 8))
        sns.heatmap(result, annot=True, fmt=".3f", cmap="coolwarm", cbar=False)
        plt.savefig(save_loc+group_name+".png")
        save_file_name= os.path.join(save_loc,group_name+".csv")
        result.to_csv(save_file_name)
        i+=1


def tt_events_statistics(tt_events_csv:str):
    tt_events=pd.read_csv(tt_events_csv)
    groups = [ "TCG1", "TCG2", "TCG3", "TCG4"]
    from scipy.stats import f_oneway
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    def perform_test(dataframe):
        # dataframe.columns = col_names
        # df_long = pd.melt(dataframe.reset_index(drop=True), var_name='Group', value_name='Value')
        # dunns_result = posthoc_dunn(df_long, group_col='Group', val_col='Value', p_adjust='bonferroni')
        # return dunns_result
        # Perform one-way ANOVA to get the F-statistic and p-value
        # anova_result = f_oneway(*[group["Value"] for name, group in df_long.groupby("Group")])
        data = dataframe.values.tolist()
        anova_result = f_oneway(data[0],data[1],data[2])
        f_statistic = anova_result.statistic
        p_value = anova_result.pvalue
        print(f"P_value = {p_value} \n f_statistic = {f_statistic}")
        # Perform Tukey's HSD test for multiple comparisons
        # tukey_result = MultiComparison(df_long['Value'], df_long['Group']).tukeyhsd()

        # Convert the Tukey's HSD test results to a DataFrame
        # tukey_df = pd.DataFrame(data=tukey_result._results_table.data[1:], columns=tukey_result._results_table.data[0])

        # Add F-statistic and p-value to the DataFrame
        # tukey_df.insert(1, "F-Statistic", f_statistic)
        # tukey_df.insert(2, "p-value", p_value)
        # return tukey_df
    # Transport events
    j = 0
    for i in range(4):
        _df = tt_events.iloc[:, j:j + 3]
        result = perform_test(_df)
        print(result)
        j +=3

def avg_std_for_bottleneck_radii():
    bottlenecks = pd.read_csv("/home/aravind/PhD_local/dean/figures/bottlenecks/"
                              "time_evolution/bottleneck_per_sim.csv")

    average = bottlenecks.mean()
    std = bottlenecks.std()

    results_df = pd.DataFrame({'Column': bottlenecks.columns, 'Average': average, 'Standard Deviation': std}).T

    results_df.to_csv("/home/aravind/PhD_local/dean/figures/bottlenecks/"
                              "time_evolution/avg_std_bottleneck.csv")


if __name__ == '__main__':
    tt_events_csv = "/home/aravind/PhD_local/dean/statistics/tt_events_HST.csv"
    tt_results = "/data/aravindramt/dean/tt/tt_0_9_5"
    unassigned_split = "/home/aravind/PhD_local/dean/figures/transport_tools/unassigned_events_sep.csv"
    groups_def = {"P1": [1, 2, 5, 7, 12, 30, 31], "P2": [3, 4, 6, 8, 11, 16, 25, 27, 41, 43, 44, 50, 58],
                  "P3": [10]}
    main_tunnel = {"P1": [1, 2, 5, 7, 12, 30, 31]}
    simulation_results = "/data/aravindramt/dean/md/simulations/"
    sim_results = "/data/aravindramt/dean/md/simulations/"
    bottleneck_csv="/home/aravind/PhD_local/dean/figures/main_images/average_bottleneck_main_figure.csv"
    details_file_loc = "/data/aravindramt/dean/tt/tt_0_9_5/data/super_clusters/details" \
                       "/initial_super_cluster_details.txt"


    # transit_time_stats(tt_results,sim_results,main_tunnel)
    # percent_events_occurance_statistics(tt_results,simulation_results,main_tunnel)
    # average_bottleneck_statistics(bottleneck_csv)
    # full_bottleneck_statistics()
    tt_events_statistics("/home/aravind/PhD_local/dean/figures/transport_tools/p1_only.csv")
    # avg_std_for_bottleneck_radii()
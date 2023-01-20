# plot bottleneck radii from transport tools output
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
# get tt results file
tt_results = pd.read_csv('/mnt/gpu/dean/tt/tt_0_9_5/statistics/4-filtered_events_statistics.txt',
                     skiprows=20, skipfooter=1, engine='python', skipinitialspace=True,
                         na_values="-")

opc_scid = [{"P1": [1], "P2": [3, 5, 6], "P3": [7, 8]}, {"P1": [1], "P2": [3, 5], "P3": [8]},
            {"P1": [1], "P2": [3, 5], "P3": [8, 24]}
    , {"P1": [1, 2, 4, 10, 11], "P2": [3, 5, 6], "P3": [7, 8]}, {"P1": [2, 4, 10, 12], "P2": [3, 9, 17, 27], "P3": []}]

tip3p_scid = [{"P1": [1], "P2": [3, 5, 6], "P3": [7, 8]}, {"P1": [1], "P2": [3, 5], "P3": [7, 8]},
              {"P1": [1, 10], "P2": [3, 5], "P3": [7, 8]}
    , {"P1": [1, 2, 4, 10, 11], "P2": [3, 5, 6], "P3": [7, 8]}, {"P1": [2, 4, 10, 12], "P2": [3, 9, 27], "P3": []}]

tip4pew_scid = [{"P1": [1], "P2": [3, 5, 6], "P3": [7, 8]}, {"P1": [1], "P2": [3, 5], "P3": [7, 8]},
                {"P1": [1], "P2": [3, 5], "P3": [8, 24]}
    , {"P1": [1, 2, 10], "P2": [6, 9], "P3": [7, 8]}, {"P1": [2, 4, 10], "P2": [3, 9], "P3": []}]
fig,ax=plt.subplots(1,3,figsize=(10,5))
for water_model in [opc_scid,tip3p_scid,tip4pew_scid]:
    tmp_df_group = pd.DataFrame()
    _bottleneck_dict = {'P1': [], 'P2': [], 'P3': []}
    for group in water_model:
        sc_ids_p1=group['P1']
        sc_ids_p2 = group['P2']
        sc_ids_p3 = group['P3']

        # Retrive the tunnel's dataframe from original result
        tmp_df_p1 = tt_results.loc[(tt_results["SC_ID"].isin(sc_ids_p1))]
        tmp_df_p2 = tt_results.loc[(tt_results["SC_ID"].isin(sc_ids_p2))]
        tmp_df_p3 = tt_results.loc[(tt_results["SC_ID"].isin(sc_ids_p3))]
        # Bottleneck

        avg_bottleneck_p1 = tmp_df_p1.Avg_BR.mean()
        avg_bottleneck_p2 = tmp_df_p2.Avg_BR.mean()
        avg_bottleneck_p3 = tmp_df_p3.Avg_BR.mean()
        _bottleneck_dict['P1'].append(avg_bottleneck_p1)
        _bottleneck_dict['P2'].append(avg_bottleneck_p2)
        _bottleneck_dict['P3'].append(avg_bottleneck_p3)
    print((_bottleneck_dict["P1"]))
    print((_bottleneck_dict["P2"]))
    print((_bottleneck_dict["P3"]))
    ax[0].bar([0,1,2,3,4],_bottleneck_dict["P1"])
    ax[1].plot(_bottleneck_dict["P2"])
    ax[2].plot(_bottleneck_dict["P3"])
plt.legend()
plt.show()



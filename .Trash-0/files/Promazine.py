# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all
#     notebook_metadata_filter: all,-language_info
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# This notebook identifies practices as part of our outlier detection who prescribed Pericyazine. The intention is that we write to them and outlying CCGs to ascertain the reasons why they use this so much compared to their peers.

# + trusted=true
#import libraries required for analysis
import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from ebmdatalab import bq, charts, maps
# -

# Data Extract Here we identify all pericyazine prescribing. We identify all prescribing to generate measures for letters.

# + trusted=true
#code from original notebook - not used as CCG codes changed 
sql = """
SELECT
  month,
  pct,
  practice,
  name,
  address1, 
  address2, 
  address4, 
  address5, 
  postcode,
  SUM(items) AS total_promazine
FROM
  ebmdatalab.hscic.normalised_prescribing_standard AS presc
INNER JOIN
  ebmdatalab.hscic.practices prac
ON
  presc.practice = prac.code
  AND (prac.setting = 4)
WHERE
  bnf_code LIKE "0402010S0%" 
###  AND
### (presc.month >= "2019-03-01" and presc.month  <= "2019-06-01")
# ## AND items > 1
GROUP BY
  month,
  pct,
  practice,
  name,
  address1, 
  address2, 
  address4, 
  address5, 
  postcode
ORDER BY
  presc.month,
  practice
"""
exportfile = os.path.join("..","data","promazine_df.csv") #set path for data cache
promazine = bq.cached_read(sql, exportfile)
promazine['month'] = promazine['month'].astype('datetime64[ns]')
promazine.head(10)
# -

# Here we identify the practices for writing to based on criteria in out methodology
#
# - prescribing in last quarter
# - at least 1 item

# + trusted=true
promazine_prescribers = promazine.loc[(promazine["month"]>= "2019-03-01") & (promazine["month"] <= "2019-06-01") & (promazine["total_promazine"] > 1)]
promazine_prescribers.head()

# + trusted=true
promazine_prescribers['practice'].nunique()

# + trusted=true
##here we identify CCGs to be written to
promazine_prescribers['pct'].unique()
# -

# we will refine this list of CCGs later

# + trusted=true
#code from old notebook - used csvs created at the time
# get data for patient list size (all patients) so we can create a measure
#sql2 = """
#SELECT month,
#pct_id AS pct,
#practice,
#sum(total_list_size) as list_size
#FROM ebmdatalab.hscic.practice_statistics
#group by 
#month,
#pct,
#practice
#order by
#month, pct
#"""
#listsize_df = bq.cached_read(sql2, csv_path='list_size.csv')
#listsize_df['month'] = listsize_df['month'].astype('datetime64[ns]')
importfile = os.path.join("..","data","list_size.csv") #set path for data cache
listsize_df = pd.read_csv(importfile)
listsize_df['month'] = listsize_df['month'].astype('datetime64[ns]')

# + trusted=true
#merge dataframes so we can do measures with deciles
promazine_per_1000 = pd.merge(left = listsize_df, right = promazine, on=['month', 'practice'], how = 'left')
promazine_per_1000['promazine_per_1000'] = 1000* (promazine_per_1000['total_promazine']/promazine_per_1000['list_size'])
promazine_per_1000['promazine_per_1000'] = promazine_per_1000['promazine_per_1000'].fillna(0)
promazine_per_1000.head(5)

# + trusted=true
exportfile = os.path.join("..","data","promazine_per_1000.csv") #set path for data cache
promazine_per_1000.to_csv(exportfile,index=False)

# + trusted=true
#create sample deciles chart
charts.deciles_chart(
        promazine_per_1000,
        period_column='month',
        column='promazine_per_1000',
        title="Items for promazine per 1000 patients (practice)",
        show_outer_percentiles=True)

#add in example https://openprescribing.net/practice/D82099/ from West Norfolk
df_subject = promazine_per_1000.loc[promazine_per_1000['practice'] == 'D82099']
plt.plot(df_subject['month'], df_subject['promazine_per_1000'], 'r--')
plt.show()

# + trusted=true
ccg_promazine = promazine_per_1000.groupby(['pct_x', 'month']).sum().reset_index()
ccg_promazine['promazine_per_1000'] = 1000* (ccg_promazine['total_promazine']/ccg_promazine['list_size'])
ccg_promazine['promazine_per_1000'] = ccg_promazine['promazine_per_1000'].fillna(0)
ccg_promazine.rename(columns = {'pct_x':'pct'}, inplace = True)
ccg_promazine = ccg_promazine.loc[(ccg_promazine["list_size"] >= 2000)]
ccg_promazine.head(5)

# + trusted=true
exportfile = os.path.join("..","data","promazine_per_1000.csv") #set path for data cache
promazine_per_1000.to_csv(exportfile,index=False)

# + [markdown]
# #ccg_pericyazine["pericyazine_per_1000"] = 1000*(ccg_pericyazine.total_pericyazine / ccg_pericyazine.list_size).fillna(0)


# + trusted=true
#create sample deciles chart
charts.deciles_chart(
        ccg_promazine,
        period_column='month',
        column='promazine_per_1000',
        title="Items for promazine per 1000 patients (CCG) ",
        show_outer_percentiles=True)

#add in example Islington is 08H
df_subject = ccg_promazine.loc[ccg_promazine['pct'] == '08H']
plt.plot(df_subject['month'], df_subject['promazine_per_1000'], 'r--')
plt.show()

# + trusted=true
#create choropeth map of cost per 1000 patients
latest_ccg = ccg_promazine.loc[(ccg_promazine["month"] == "2017-08-01")]
plt.figure(figsize=(12, 7))
plt = maps.ccg_map(
    latest_ccg, 
    title="Items for promazine per 1000 patients (CCG) ", 
    map_year = '2018',
    column='promazine_per_1000', separate_london=False,
    plot_options={'cmap': 'coolwarm'}
    )
plt.show()

# + trusted=true
## here we look at CCGs that are outliers.
latest_ccg.sort_values("total_promazine", ascending=False).head(250)

# + trusted=true
exportfile = os.path.join("..","data","latest_ccg_promazine.csv") #set path for data cache
latest_ccg.to_csv(exportfile,index=False)
# -

print(exportfile)

# + trusted=true
import glob
import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path


def ccg_map_bespoke(
    df,
    title="",
    column=None,
    region=None,
    separate_region=False,
    london_layout="horizontal",
    cartogram=False,
    subplot_spec=None,
    show_legend=True,
    map_year=None,
    plot_options=None,
):
    """Draw a CCG map with London separated out
    """
    # Because this uses subplots to arrange London and England,
    # the only way to arrange nested subplots is via a subplotspec
    assert column, "You must specify a column name to plot"
    df = df.copy()
    # input df must have 'pct' column, plus others as specified
    data_dir = os.path.join("..","data")

    # Add names and if it's London. Note the names in ccg_for_map must
    # match names in the CCG geojson, as that doesn't include codes at
    # the momemt
    map_name = os.path.join(data_dir, "ccg_for_map.csv")
    names = pd.read_csv(map_name)

    # Check we know about all the codes in the input data
    diff = np.setdiff1d(df["pct"], names["code"])
    if len(diff) > 0:
        raise BaseException(
            "Data contains CCG codes we don't know about: {}".format(diff)
        )

    df = df.merge(names[["code", "name", "region"]], left_on="pct", right_on="code")
    df = df.set_index("name")

    # Load map data
    cartogram_suffix = ""
    if cartogram:
        cartogram_suffix = "_cartogram"
    if map_year:
        map_file = os.path.join(data_dir,"ccgs{}_{}.json").format(cartogram_suffix, map_year)
    else:
        map_file = sorted(
            glob.glob(str(os.path.join(data_dir,"ccgs{}_2*.json").format(cartogram_suffix)))
        )[-1]
    ccgs = gpd.read_file(map_file)
    # Normalise names to match `ccg_fo_map` format (above)
    ccgs["name"] = ccgs["name"].str.upper()
    ccgs = ccgs.set_index("name")
    # Remove ones without geometry - these are (usually) federations
    # rather than individual CCGs
    ccgs = ccgs[~ccgs["geometry"].isnull()]

    # Check we can map all the CCGs named in the input data
    diff = np.setdiff1d(df.index, ccgs.index)
    if len(diff) > 0:
        raise BaseException("Data contains CCG names we can't map: {}".format(diff))

    # Join map with data
    gdf = ccgs.join(df, rsuffix="_orig")

    # Split into london and rest of England
    gdf_london = gdf[gdf["region"] == region]
    gdf_roe = gdf

    # set common value limits for colour scale
    default_plot_options = {
        'vmin': df[column].min(),
        'vmax': df[column].max(),
        'edgecolor': "black",
        'linewidth': 0.1,
        'cmap': "OrRd",
    }

    if plot_options is None:
        plot_options = {}

    for k, v in default_plot_options.items():
        if k not in plot_options:
            plot_options[k] = v

    def plot(gdf, ax, title="", legend=True):
        gdf.plot(
            ax=ax,
            column=column,
            legend=legend,
            **plot_options
        )
        ax.set_aspect(1.63)
        if title:
            ax.set_title(title, size=14)
        ax.axis("off")

    fig = plt.gcf()
    if not subplot_spec:
        subplot_spec = gridspec.GridSpec(1, 1)[0]
    if separate_region:
        if london_layout == "horizontal":
            gs = gridspec.GridSpecFromSubplotSpec(
                nrows=1, ncols=2, width_ratios=[1, 2], subplot_spec=subplot_spec
            )
            ldn_ax = fig.add_subplot(gs[0, 0])
            roe_ax = fig.add_subplot(gs[0, 1])
        else:
            gs = gridspec.GridSpecFromSubplotSpec(
                nrows=2, ncols=1, height_ratios=[2, 1], subplot_spec=subplot_spec
            )
            roe_ax = fig.add_subplot(gs[0, 0])
            ldn_ax = fig.add_subplot(gs[1, 0])

        plot(
            gdf_roe,
            roe_ax,
            title="England:\n{}".format(title),
            legend=show_legend,
        )
        plot(gdf_london, ldn_ax, title=":\n{}".format(title), legend=False)
    else:
        ax = plt.subplot(subplot_spec)
        plot(gdf, ax, title=title, legend=show_legend)
    return plt


# + trusted=true
plt = ccg_map_bespoke(
    latest_ccg, 
    title="Items for promazine per 1000 patients (CCG) ", 
    map_year = '2018',
    column='promazine_per_1000', region='North West', separate_region=True,
    plot_options={'cmap': 'coolwarm'}
    )

# + trusted=true



# + trusted=true
    names = os.path.join("..","data","ccg_for_map_eoe.csv")


# + trusted=true
print(names)
# -



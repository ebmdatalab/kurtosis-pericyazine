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
#   vscode:
#     interpreter:
#       hash: de1343822d6e7d7aeea8796be9d48304b0fa3610166e8740495ec86b33c71a9e
# ---

# # Kurtosis and similar methods for identifying outliers in prescribing

# This notebook looks at each chemical prescribed in England (using national prescribing data), assigns it to its class (BNF subparagraph) and, for each Clinical Commissioning Group, calculates the chemical-class proportion. 
#
# Chemicals are then ranked by (a) their kurtosis and (b) a ratio between inter-centile differences, in order to identify those with extreme distributions, i.e. outliers. 

# + trusted=true
import pandas as pd
import numpy as np
from ebmdatalab import bq
import os

q = '''SELECT p.* FROM ebmdatalab.outlier_detection.chem_by_subpara_by_ccg_juntoaug17_limitsubpara p
-- exclude non-standard CCG codes:
INNER JOIN ebmdatalab.hscic.ccgs c ON p.pct = c.code AND c.org_type = "CCG"
'''
# see chemical_by_subpara.sql for data source production

df1 = bq.cached_read(q, csv_path=os.path.join("..","data","df1.csv"), use_cache=True)
# rows: pct, chemical, subpara, num, denom, ratio (num and denom are items not quantity)

q2 = '''SELECT DISTINCT chemical, chemical_code from ebmdatalab.hscic.bnf'''
chem = bq.cached_read(q2, csv_path=os.path.join("..","data","chemical.csv"), use_cache=True)

q3 = '''SELECT DISTINCT subpara, subpara_code from ebmdatalab.hscic.bnf'''
subp = bq.cached_read(q3, csv_path=os.path.join("..","data","subpara.csv"), use_cache=True)

df1.head()

num_CCGs = df1.pct.nunique()

# -

# ## Clean the data: sort out null denominators

# + trusted=true
# need to flag where ccgs have not prescribed any items of the denominator in order to clean the data. 

# Step 1: amend the datafrome to include a line for every CCG and every chemical and subparagraph,, even if not prescribed.

# list all subpara-chemical combinations 
a = df1[["subpara", "chemical"]].drop_duplicates()

#list all ccgs
b = df1[["pct"]].drop_duplicates()

# cross join to make table of all CCGs and all subpara combinations 
a['tmp'] = 1
b['tmp'] = 1
c = b.merge(a, on="tmp").drop('tmp', axis=1) # 237,636 rows

# join to data - this will list every possible chemical against every CCG
data = c.merge(df1, how="left", on=["pct","subpara","chemical"])  # 237,636 rows
data


# Step 2: identify those with zero subparas
# subpara totals by ccg
subpara = df1[["pct","subpara","denom"]].groupby(["subpara","pct"]).max().reset_index() # 42,917 rows

#list all possible subparagraphs and all ccgs
a2 = df1[["subpara"]].drop_duplicates()
a2['tmp'] = 1

# cross join to CCGs to make table of all CCGs and all subpara combinations 
c2 = b.merge(a2, on="tmp").drop('tmp', axis=1) # 56,097 rows

# join to subpara data by ccg to identify subparas prescribed by each ccg.  
d = c2.merge(subpara,how="left", on=["subpara","pct"])

# for subparas never prescribed, replace NAs with zeros so that there is data present to indicate this
d = d.fillna(0)

# join back to original dataset
d2 = d.merge(data, how="left", on=["subpara","pct"], suffixes=("_subpara",""))
# check how many have zero denominators:
# data.loc[(data["denom_subpara"]==0)] # 19,665 rows 

# exclude combinations where denominators are zero THEN replace NAs with 0:
data2 = d2.loc[(d2["denom_subpara"]!=0)]
data2 = data2.fillna(0)
data2.head()

# print( f"Number of CCGs in data2: {data2.pct.nunique()}")
# -

# ### Filter out low numbers (chemical and subpara)

# + trusted=true
# total prescribing for each chemical
# sum numerators to find total volume for each chemical across all ccgs
num = pd.DataFrame(df1["num"].groupby(df1["chemical"]).sum()).reset_index()

# calculate centile for numerator for each ccg
num["num centile"] = pd.qcut(num["num"], 10, labels=np.arange(1,11,1))

# total prescribing for each paragraph
d3 = d2[["pct","subpara","denom_subpara"]].drop_duplicates()
d3 = d3.groupby("subpara").sum().sort_values(by="denom_subpara")
d3["denom centile"] = pd.qcut(d3["denom_subpara"], 10, labels=np.arange(1,11,1))
d3 = d3.reset_index()


# merge with data table
data3 = data2.merge(num, how="inner", on="chemical",suffixes=("","_total"))
data3 = data3.merge(d3, how="inner", on="subpara",suffixes=("","_total"))

num_chemicals_stage0 = data3['chemical'].nunique()

# filter out CCGs with lowest denominator values
national_count_minimum = 1000
data3 = data3.loc[(data3["denom centile"]>2) & (data3["num_total"]>national_count_minimum)]

num_chemicals_stage1 = data3['chemical'].nunique()
num_chemicals_lost_01 = num_chemicals_stage0 - num_chemicals_stage1

# num of CCGs prescribing each chemical
count_ = pd.DataFrame(data3.loc[data3["num"]>0].groupby("chemical")["pct"].nunique()).reset_index()
count_ = count_.rename(columns={"pct":"count2"})

data3 = data3.merge(count_, how="inner", on="chemical")

data3.head()

# -

# ## Calculate key stats
# ### Median, Range, SD, Kurtosis and Skew

# + trusted=true
#select columns of interest and get key stats
df2 = pd.DataFrame(data3.groupby(["chemical","subpara","num_total","num centile","denom_subpara_total","denom centile","count2"])["ratio"].describe())
#df2 = df2.unstack()
#df2.columns = df2.columns.droplevel()

df3 = df2.reset_index()
df3["range"] = df3["max"] - df3["min"]
df3 = df3[["chemical","subpara","num_total","num centile","denom_subpara_total","denom centile","count","count2","50%","min","max","range","std"]].rename(columns={"50%":"median"})

num_chemicals_stage2 = df3['chemical'].nunique()

# filter out chemicals in paragraphs only prescribed by few CCGs
df3 = df3.loc[df3["count"]>=50]

num_chemicals_stage3 = df3['chemical'].nunique()
num_chemicals_lost_23 = num_chemicals_stage2 - num_chemicals_stage3


# reshape data to put CCGs in columns
df5 = data3.pivot(index="chemical",columns='pct', values='ratio')

#calculate kurtosis and skew for each chemical
import scipy.stats as stats
k = pd.Series(stats.kurtosis(df5, axis=1,nan_policy="omit"),name="kurtosis")
sk =  pd.Series(stats.skew(df5, axis=1,nan_policy="omit"),name="skew")

# compile all stats together
result = pd.concat([df3, k, sk], axis=1).sort_values(by="kurtosis",ascending=False)
result = result[["chemical","subpara", "count","count2","num_total","num centile","median","min","max","range","std","kurtosis","skew"]].round(2)


# Lookup chemical and subparagraph names
df4 = result.merge(chem, how="left", left_on = "chemical",right_on="chemical_code",suffixes=(""," name"))
df4 = df4.merge(subp, how="left", left_on = "subpara",right_on="subpara_code",suffixes=(""," name"))
df4 = df4[["chemical","chemical name","subpara","subpara name","num_total","num centile", "count","count2","median","min","max","range", "std","kurtosis","skew"]].round(2)
df4["subpara"] = df4["subpara"].fillna(0).astype(int).astype(str)

df4.head()

# + trusted=true
### df5 is a chemical X CCG - raw data

print( df5.head )
print( df5.columns )
print( df5.shape )


# + trusted=true
### df4 contains df5 data, processed and summarised, including kurtosis and skew

print( df4.head )
print( df4.columns )
print( df4.shape )

# -

# ## Ranking Chemicals by Range, Kurtosis, Skew and SD
# Those with high SD tend to be chemicals where there is general disagreement, so outliers are not *that* unusual

# + scrolled=true trusted=true
# limit to those with positive skew, i.e. most CCGs prescribe few and those prescribing more are ouliers,
# and range at least 10%:

num_chemicals_stage4 = df4['chemical'].nunique()

dfp = df4.loc[(df4["skew"]>=0) & (df4["range"]>0.1)]

num_chemicals_stage5 = dfp['chemical'].nunique()
num_chemicals_lost_45 = num_chemicals_stage4 - num_chemicals_stage5


# sort by range
r1 = dfp.sort_values(by=["range","kurtosis"],ascending=False)
# create a ranking
r1["R"]  = r1["range"].rank(ascending=False, method="min")

# sort by kurtosis 
r2 = dfp.sort_values(by=["kurtosis"],ascending=False)
r2["K"] = r2["kurtosis"].rank(ascending=False, method="min")

# sort by skew
r3 = dfp.sort_values(by=["skew"],ascending=False)
r3["Sk"] = r3["skew"].rank(ascending=False, method="min")

r1 = dfp.copy()
# create a ranking
r1["R"] = r1["range"].rank(ascending=False, method="min")
r1["K"] = r1["kurtosis"].rank(ascending=False, method="min")
r1["Sk"] = r1["skew"].rank(ascending=False, method="min")


r1.sort_values(by="K").head(20)

# + trusted=true
### dfp filters df4 for only those chemicals showing a positive skew and a range of >10%

print( dfp.head )
print( dfp.columns )
print( dfp.shape )

# -

# ### Add hyperlinks to maps

# + trusted=true
# assign overall scores based on all three rankings

r2 = r1.copy()
columns = ["R","K","Sk"]
r2['score'] = 0
for col in columns:
    r2.loc[r2[col]<=10, 'score'] = r2['score']+10 # if in top 10, score 10
    r2.loc[(r2[col]>10)&(r2[col]<=50), 'score'] = r2['score']+5 # if in top 10, score 10
rc = r2.sort_values(by=["score","kurtosis"],ascending=False)
rc.head()
# -

# note: links go to up-to-date maps and can only show chemical/para not subpara. 

# + trusted=true
# NBVAL_IGNORE_OUTPUT
# ^this is a magic comment to help tests pass

links = rc.copy().head(50)
links["str"] = links["subpara"].str[0:2].map(int).map(str)+"."+links["subpara"].str[2:4].map(int).map(str)+"."+links["subpara"].str[4:6].map(int).map(str)
links["link"] = "https://openprescribing.net/analyse/#org=CCG&numIds="+links["chemical"]+"&denomIds="+links["str"]+"&selectedTab=map"
links = links.drop("str",axis=1)

def make_clickable(val):
    return '<a href="{}">{}</a>'.format(val,val)

links = pd.DataFrame(links).style.format(make_clickable, subset=['link'])
links
# -

# ### Histograms

# + trusted=true
# plot top 50 only
## exclude practice IDs etc and select only lists of figures to plot:
dfh = rc[["chemical","chemical name","subpara name","count","count2","score","kurtosis","num_total"]].head(50).merge(data3[["chemical","ratio"]], on="chemical").drop("chemical",axis=1)#.sort_values(by="chemical name")

# create list of titles for charts - combine chem name with para name and no of CCGs prescribing each 
titles = dfh.groupby(["chemical name","subpara name","kurtosis","count","count2","num_total"]).count().reset_index().rename(columns={"ratio":"n"})
titles["title"] = titles["chemical name"]+", K="+titles["kurtosis"].round(1).map(str)+"\n"+titles["count2"].map(int).map(str)+" CCGs prescribed "+titles["num_total"].map(int).map(str)+" items\n"+titles["subpara name"].str[0:35]+" ("+titles["count"].map(int).map(str)+" CCGs)"

titles = dfh[["chemical name"]].drop_duplicates().merge(titles,on="chemical name")
titles = titles["title"]



import matplotlib.pyplot as plt
import seaborn as sns
from textwrap import wrap # this will interpret "\n" as newline

## use facetgrid to create a grid of charts. 
# "chemical name" column is used to split data across charts. 
#### note it's also possible to set hue="columnname" to show multiple data groupings on one chart split by colour. 
g = sns.FacetGrid(dfh, col="chemical name",col_wrap=3,sharey=True,sharex=False,height=5)

## define histograms to plot:
g.map(sns.distplot, "ratio", kde=False, bins=20, hist_kws={"color": '#81c5e2', "alpha": .6}) 
g.set_ylabels("Frequency (CCGs)", size=15, alpha=0.5)
g.set_xlabels("Proportion of Subpara", size=15, alpha=0.5)

## Set the ticklabel size and color:
#loop through each plot
for ax,title in zip(g.axes.flat,titles):
    ax.tick_params(labelsize=10,labelcolor="black")
    
## Set the size of titles:
g.set_titles(size = 22)

## Set the titles
for ax, title in zip(g.axes.flat, titles):
    ax.set_title(title)

## adjust spacing of sub plots:
plt.subplots_adjust(left=0.15, top=0.8, hspace = 0.3)

plt.show()
# -

# ## Ranking chemicals by percentile differences / IQR etc.
# First Calculate various percentiles

# + trusted=true
dftest = df5.transpose()
q = [0.03,0.05,0.25,0.5,0.75,0.95,0.97]

# calculate modal value for each chemical
# note where there are multiple modal values, all are given so we take max here
mo = pd.DataFrame(dftest.mode(numeric_only=True).max()).reset_index()
mo = mo.rename(columns={0:"mode"})

smy0 = dftest.describe(percentiles=q).drop(["count","mean","std"]).transpose().reset_index()

#calculate IQR
smy0["IQR"] = smy0["75%"]-smy0["25%"] 
smy0 = smy0.merge(mo, on="chemical", how="left")

smy0.head()

# -

# ### Chemicals with largest ratio of 95-97th percentiles to 50-95th percentiles and with mode = 0
# A large difference between the top-prescribing CCGs and the mid-high prescribing CCGs will indicate that there are several CCGs prescribing well above average levels.

# + trusted=true

smy0[['mode','50%']].plot.kde()
plt.axvline(x=0,c="grey",ls='--')
smy0[['mode','50%']].hist()

smy0[['mode']].min()
smy0[['50%']].min()


# + trusted=true
smy2 = smy0.copy()
smy2["95-97"] = smy2["97%"]-smy2["95%"]
smy2["50-95"] = smy2["95%"]-smy2["50%"]
smy2["ratio2"] = smy2["95-97"]/smy2["50-95"]

num_chemicals_stage6 = smy2['chemical'].nunique()

# limit to those where mode is zero
num_mode = sum(smy2["mode"]==0)
num_median = sum(smy2["50%"]<0.1)
print( f"Using the mode: {num_mode}" )
print( f"Using the median: {num_median}" )

# smy2 = smy2.loc[smy2["mode"]==0]
smy2 = smy2.loc[smy2["50%"]<0.1]

num_chemicals_stage7 = smy2['chemical'].nunique()
num_chemicals_lost_67 = num_chemicals_stage6 - num_chemicals_stage7


smy2 = smy2.reset_index()
smy2 = smy2[["chemical","ratio2","25%"]].merge(dfp, on="chemical")
smy2["M2"]  = smy2["ratio2"].rank(ascending=False, method="min")

smy2.sort_values(by="ratio2",ascending=False).head(10)



num_chemicals_stage8 = smy2['chemical'].nunique()
num_chemicals_lost_78 = num_chemicals_stage7 - num_chemicals_stage8


# + trusted=true
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D  # for legend handle

smy2["r_rank"]  = smy2["ratio2"].rank(ascending=False, method="min")
smy2["k_rank"]  = smy2["kurtosis"].rank(ascending=False, method="min")
num_chemicals_plotted = smy2['chemical'].nunique()

n_to_label = 5

data_to_plot = smy2[['ratio2','kurtosis','chemical name','r_rank','k_rank']].rename(columns={"chemical name":"label", "ratio2":"ratio"})
data_to_plot.set_index('label', inplace=True)


## Which chemicals to highlight?
data_to_plot['group'] = f"Rank > {n_to_label}"
data_to_plot.loc[(data_to_plot.r_rank<=n_to_label)|(data_to_plot.k_rank<=n_to_label), 'group'] = 'Rank <= 5'
colour_map = {f"Rank <= {n_to_label}":'orange', f"Rank > {n_to_label}":'grey'}
alpha_map = {f"Rank <= {n_to_label}":1, f"Rank > {n_to_label}":0.4}

fig, ax = plt.subplots()

ax.scatter(data_to_plot['ratio'], data_to_plot['kurtosis'], c=data_to_plot['group'].map(colour_map), alpha=data_to_plot['group'].map(alpha_map) )

for k, v in data_to_plot.iterrows():
    if v.group == f"Rank <= {n_to_label}":
        ax.annotate(k, (v['ratio'], v['kurtosis']), textcoords="offset points", xytext=(5,-15), ha='left', rotation=30)

handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=v, label=k, markersize=8) for k, v in colour_map.items()]
ax.legend(title='', handles=handles, loc='lower right')
plt.xlabel("high:mid centile ratio", weight="bold")
plt.ylabel("Kurtosis", weight="bold")
plt.suptitle(f"Outlier metrics for all candidate chemicals (n={num_chemicals_plotted})", weight="bold", x=0.125, y=0.99, horizontalalignment="left")
plt.title(f"high:mid centile ratio and kurtosis", x=0, y=1.04, horizontalalignment='left' )
plt.gcf().set_size_inches(8, 6)

plt.savefig('ratio_kurtosis_plot.png', dpi=200, bbox_inches='tight')

plt.show()




# + trusted=true

print( smy2.loc[smy2['chemical name'] == 'Pericyazine', ['chemical name', 'r_rank', 'k_rank']] )
print( smy2.loc[smy2['chemical name'] == 'Promazine hydrochloride', ['chemical name', 'r_rank', 'k_rank']] )


# + trusted=true
x = data_to_plot['ratio']
y = data_to_plot['kurtosis']
c = data_to_plot['group'].map(colour_map)
a = alpha=data_to_plot['group'].map(alpha_map)

def scatter_hist(x, y, c, a, ax, ax_histx, ax_histy, x_bins, y_bins):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y, color=c, alpha=a)

    ax_histx.hist(x, bins=x_bins)
    ax_histy.hist(y, bins=y_bins, orientation='horizontal')
    

# Start with a square Figure.
fig = plt.figure(figsize=(6, 6))
# Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
# the size of the marginal axes and the main axes in both directions.
# Also adjust the subplot parameters for a square plot.
gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)
# Create the Axes.
ax = fig.add_subplot(gs[1, 0])
plt.xlabel("high:mid centile ratio", weight="bold")
plt.ylabel("Kurtosis", weight="bold")
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
plt.ylabel("Fequency", weight="bold")
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
plt.xlabel("Fequency", weight="bold")

# Draw the scatter plot and marginals.
scatter_hist(x, y,
             c, a,
             ax, ax_histx, ax_histy,
             x_bins = np.arange(0,2,0.05),
             y_bins = np.arange(0,130,5) )


for k, v in data_to_plot.iterrows():
    if v.group == f"Rank <= {n_to_label}":
        if ( k in ['Tolbutamide', 'Cefadroxil'] ):
            ax.annotate(k, (v['ratio'], v['kurtosis']), textcoords="offset points", xytext=(0,15), ha='center' )
        elif ( k in ['Mefloquine hydrochloride'] ):
            ax.annotate(k, (v['ratio'], v['kurtosis']), textcoords="offset points", xytext=(0,15), ha='left' )
        elif ( k in ['Pericyazine'] ):
            ax.annotate(k, (v['ratio'], v['kurtosis']), textcoords="offset points", xytext=(0,-20), ha='center', weight="bold" )
        elif ( k in ['Levofloxacin'] ):
            ax.annotate(k, (v['ratio'], v['kurtosis']), textcoords="offset points", xytext=(0,-20), ha='center' )
        elif ( k in ['Other vitamin B compound preparations'] ):
            ax.annotate(k, (v['ratio'], v['kurtosis']), textcoords="offset points", xytext=(-20,-20), ha='left' )
        elif ( k in ['Isotretinoin','Balsalazide sodium'] ):
            ax.annotate(k, (v['ratio'], v['kurtosis']), textcoords="offset points", xytext=(10,0), ha='left' )
        elif ( k in ['Promazine hydrochloride'] ):
            ax.annotate(k, (v['ratio'], v['kurtosis']), textcoords="offset points", xytext=(30,-20), ha='right', weight="bold")
        else:
            ax.annotate(k, (v['ratio'], v['kurtosis']), textcoords="offset points", xytext=(20,-20), ha='right', rotation=30)

handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=v, label=k, markersize=8) for k, v in colour_map.items()]
ax.legend(title='', handles=handles, loc=[1.029,1.03])
# plt.xlabel("high:mid centile ratio", weight="bold")
# plt.ylabel("Kurtosis", weight="bold")
plt.suptitle(f"Outlier metrics for all candidate chemicals (N={num_chemicals_plotted})", weight="bold", x=0.1, y=0.95, horizontalalignment="left")
# plt.title(f"high:mid centile ratio and kurtosis", x=0, y=1.04, horizontalalignment='left' )
plt.gcf().set_size_inches(8, 6)

plt.savefig('ratio_kurtosis_plot_with_histograms.png', dpi=300, bbox_inches='tight')

plt.show()
# -

# ### Histograms for top chemicals by percentile ratio (95-97th:50-95th) and with mode = 0

# + trusted=true
## exclude practice IDs etc and select only lists of figures to plot:
dfh = smy2.sort_values(by="ratio2",ascending=False).head(25)
dfh = dfh[["chemical","chemical name","subpara name","kurtosis","ratio2"]].merge(data3[["chemical","ratio"]], on="chemical")#.sort_values(by="chemical name")

# create list of titles for charts - combine chem name with para name and no of CCGs prescribing each 
titles = dfh.groupby(["chemical","chemical name","subpara name","kurtosis","ratio2"]).count().reset_index().rename(columns={"ratio":"n"})
titles["title"] = titles["chemical"]+"\n"+ titles["chemical name"]+" K="+titles["kurtosis"].map(str)+"\n ("+titles["subpara name"].str[0:35]+", n="+titles["n"].map(str) +")"
titles = dfh[["chemical name"]].drop_duplicates().merge(titles,on="chemical name")
titles = titles["title"]



import seaborn as sns
from textwrap import wrap # this will interpret "\n" as newline

## use facetgrid to create a grid of charts. 
# "chemical name" column is used to split data across charts. 
#### note it's also possible to set hue="columnname" to show multiple data groupings on one chart split by colour. 
g = sns.FacetGrid(dfh, col="chemical name",col_wrap=3,sharey=True,sharex=False,height=5)

## define histograms to plot:
g.map(sns.distplot, "ratio", kde=False, bins=20, hist_kws={"color": '#81c5e2', "alpha": .6}) 
g.set_ylabels("Frequency (CCGs)", size=15, alpha=0.5)
g.set_xlabels("Proportion of Subpara", size=15, alpha=0.5)

## Set the ticklabel size and color:
#loop through each plot
for ax,title in zip(g.axes.flat,titles):
    ax.tick_params(labelsize=10,labelcolor="black")
    
## Set the size of titles:
g.set_titles(size = 22)

## Set the titles
for ax, title in zip(g.axes.flat, titles):
    ax.set_title(title)

## adjust spacing of sub plots:
plt.subplots_adjust(left=0.15, top=0.8, hspace = 0.3)

plt.show()

# + trusted=true
print(f"--- Chemical filter stage 0 -> 1 ---")
print(f"Starting with {num_chemicals_stage0} chemicals")
print( "Filtering by decile (>2)")
print(f"{num_chemicals_lost_01} chemicals are removed from the bottom two deciles")
print(f"We are left with {num_chemicals_stage1} chemicals")
print(f"--- Chemical filter stage 2 -> 3 ---")
print(f"Starting with {num_chemicals_stage2} chemicals")
print( "Filtering by CCG numbers (>=50)")
print(f"{num_chemicals_lost_23} chemicals are removed due to low CCG numbers")
print(f"We are left with {num_chemicals_stage3} chemicals")
print(f"--- Chemical filter stage 4 -> 5 ---")
print(f"Starting with {num_chemicals_stage4} chemicals")
print( "Filtering by positive skew and range > 10%")
print(f"{num_chemicals_lost_45} chemicals are removed due to skew/range filtering")
print(f"We are left with {num_chemicals_stage5} chemicals")
print(f"--- Chemical filter stage 6 -> 7 ---")
print(f"Starting with {num_chemicals_stage6} chemicals")
print( "Filtering by modal proportion == 0")
print(f"{num_chemicals_lost_67} chemicals are removed due to modal proportion constraints")
print(f"We are left with {num_chemicals_stage7} chemicals")
print(f"--- Chemical filter stage 7 -> 8 ---")
print(f"Starting with {num_chemicals_stage7} chemicals")
print( "Merging with dfp")
print(f"{num_chemicals_lost_78} chemicals are removed after dfp merge")
print(f"We are left with {num_chemicals_stage8} chemicals")


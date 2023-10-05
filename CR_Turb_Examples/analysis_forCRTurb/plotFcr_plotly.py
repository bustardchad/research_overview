import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
import time
import seaborn as sns


def add_legend(plot_args_lst, legend_loc = 'best', default_plot_args = {'color' : 'black'},
                 **kwargs):
      """Adds another legend to plot.
      Keyword Arguments:
      plot_args_lst      -- List of plotting arguments to show.
      legend_loc         -- Location of new legend (default 'best')
      default_plot_args  -- Arguments will be used for every item in new legend.
      kwargs             -- Will be passed to `plt.legend` of new legend.
      Example:
              add_legend([{'ls' : '-', 'label' : r'$\rho > \rho_{cl}/3$'}, {'ls' : '--', 'label' : r'$T   < 2 T_{cl}$'}],
                                  default_plot_args = {'c' : 'k'},
                                  title = r'')

      Will add a legend with two different lines (both black).
      """
      ax = plt.gca()
      leg = ax.get_legend()
      linelst = []
      for cargs in plot_args_lst:
          for k, v in default_plot_args.items():
              if k not in cargs:
                  cargs[k] = v
          l, = plt.plot(np.nan, **cargs)
          linelst.append(l)
      o = kwargs.copy()
      if 'loc' not in o:
          o['loc'] = legend_loc
      if 'fontsize' not in o:
          o['fontsize'] = 12
      if legend_loc == 'above':
          o['loc'] = 'lower center'
          o['bbox_to_anchor'] = (0.5, 1.01)
      plt.legend(handles = linelst, **o)
      if leg is not None:
          ax.add_artist(leg) # Add old legend



#####################################################3
# Simulation details:
#   CR transport is diffusion only
#   Varying pc/pg, keeping beta ~ 1 constant
#   Run on 256^3 box, driving turbulence at k = 2 mode

# Isothermal equation of state
pcpg_iso_Mach05 = [1e-4,5e-2,0.2, 0.3, 1, 1.5, 10, 100]
pcpg_iso_Mach035 = [1.5,10,100]

fcr_iso_Mach05 = [9.6e-4,0.14, 0.58, 0.65, 0.8, 0.86, 0.89, 0.86]
fcr_iso_Mach035 = [0.83, 0.91, 0.85]


# Same but adiabatic equation of state instead of isothermal
pcpg_ad_Mach05 = [0.025,0.06,0.1,1/5,1/2,1,10,100]
fcr_ad_Mach05 = [0.12,0.19,0.23,0.4,0.56, 0.68, 0.84, 0.79]
fg_ad_Mach05 = [0.78,0.7,0.62,0.47,0.3,0.21, 0.07, 0.14]


# adiabatic eos, now with no CR diffusion (only advective transport)
pcpg_advect_Mach05 = [1e-2,1e-1,1,10,100]
fcr_advect_Mach05 = [0.0003,0.0003,0.096,0.52,0.73]
fg_advect_Mach05 = [0.969,0.969,0.74,0.35,0.13]

# expectation given diffusion-only CR reacceleration rates (see Bustard and Oh 2022b)
#pcpg_arr = np.arange(1e-4,2.0,0.0001)
pcpg_arr = np.arange(1.e-4,1.e2,0.0001)
fcr_expect = (2./3.)*(5e6/(1e7*np.sqrt(1.0 + pcpg_arr)))*(pcpg_arr*1e7**2.0)/(5e6**2.0)
fcr_expect_referee = fcr_expect/(1.0 + fcr_expect)


fill_style='full'
ms = 12.0
plt.semilogx(pcpg_ad_Mach05,fcr_ad_Mach05,'ko-',label=r"$f_{CR}$",markersize=ms,fillstyle=fill_style)
plt.semilogx(pcpg_ad_Mach05,fg_ad_Mach05,'ro-',label=r"$f_{th}$",markersize=ms,fillstyle=fill_style)
plt.semilogx(pcpg_advect_Mach05,fcr_advect_Mach05,'kD-',markersize=ms,fillstyle=fill_style)
plt.semilogx(pcpg_advect_Mach05,fg_advect_Mach05,'rD-',markersize=ms,fillstyle=fill_style)
plt.semilogx(pcpg_arr,fcr_expect,'b-',label=r"$\frac{2}{3} \frac{M_{ph} P_{CR}}{\rho v^{2}}$",markersize=ms,fillstyle=fill_style)
plt.ylim(8e-3,1)
plt.xlim(1e-2,125)
plt.xlabel(r"$P_{CR}/P_{g}$",fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title("Diffusion Only",fontsize=22)
plt.text(0.015,0.88,r"$\beta \sim 1$",fontsize=16,bbox=dict(edgecolor='black', alpha=0.1))
plt.legend(prop={'size': 12},ncol=2,bbox_to_anchor=(0.4, -0.12),title=r'$\bf{Color}$',frameon=False)
# add another legend to distinguish diffusion vs no diffusion cases
add_legend([{'marker' : 'o', 'label' : r'$\kappa \sim 0.15 L_{0}v_{ph}$'}, {'marker' : 'd', 'label' : r'$\kappa = 0$'}],
         default_plot_args = {'ls' : '', 'c' : 'k', 'markersize' : '12','fillstyle' : 'full'}, ncol=1,bbox_to_anchor=(1.02, -0.12),frameon=False,
                                   title = r'$\bf{Symbol}$')
plt.tight_layout()
plt.savefig("fcr_vs_pcpg_adiabatic.pdf")
plt.close()




# that plot is ugly...let's make it a stack plot instead
y_diff = np.vstack([fcr_ad_Mach05,fg_ad_Mach05])
y_nodiff = np.vstack([fcr_advect_Mach05,fg_advect_Mach05])

labels = [r"f$_{\rm CR}$", r"f$_{\rm th}$"]
color_map=["blue","orange"]

fig, axs = plt.subplots(2,1,sharex=True,sharey=True,figsize=(5,6))
# top plot is with diffusive CR transport
axs[0].stackplot(pcpg_ad_Mach05,y_diff,labels=labels,colors=color_map,alpha=0.5,lw=1.5, edgecolor='black')
axs[0].plot(pcpg_arr,fcr_expect,'k--',label=r"$\frac{2}{3} \frac{\mathcal{M}_{ph} P_{\rm CR}}{\rho v^{2}}$",linewidth=2)
axs[0].plot(pcpg_arr,fcr_expect_referee,'k-.',label=r"$\frac{\frac{2}{3} \frac{\mathcal{M}_{ph} P_{\rm CR}}{\rho v^{2}}}{(1 + \frac{2}{3} \frac{\mathcal{M}_{ph} P_{\rm CR}}{\rho v^{2}})}$",linewidth=2)
plt.xscale('log')
axs[0].set_ylabel(r"$f_{i}$",fontsize=15)
#axs[0].legend(loc='lower right',title=r"$\kappa_{||} \sim 0.15 L_0 c_s$")
#axs[0].legend(loc='upper left')
#plt.ylabel(r"$\dot{E}/\epsilon$",fontsize=16)
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center',ncol=4,prop={'size': 10.0})
# bottom plot is no diffusion
axs[1].stackplot(pcpg_advect_Mach05,y_nodiff,labels=labels,colors=color_map,alpha = 0.5,lw=1.5, edgecolor='black')
axs[1].set_xlim(min(pcpg_ad_Mach05),max(pcpg_ad_Mach05))
axs[1].set_ylim(0,1)

props = dict(boxstyle='round', facecolor='white', alpha=0.0)
axs[0].text(0.57, 0.6, "CRs back-react on" +"\n" + "flow, " + r"$f_{\rm CR}$ saturates", transform=axs[0].transAxes, fontsize=12,verticalalignment='top', bbox=props)
#mline_label = r"$f_{\rm CR}$ follows" + "\n analytic expectation"
#axs[0].text(0.44, 0.24, mline_label, transform=axs[0].transAxes, fontsize=12,
            #verticalalignment='top',horizontalalignment='center', bbox=props)
mline_label2 = r"$f_{\rm CR} > 0$ due to" + "\n numerical dissipation"
axs[1].text(0.75, 0.28, mline_label2, transform=axs[1].transAxes, fontsize=12,
            verticalalignment='top',horizontalalignment='center', bbox=props)


plt.xscale('log')
plt.xlabel(r"$P_{CR}/P_{g}$",fontsize=15)
axs[1].set_ylabel(r"$f_{i}$",fontsize=15)
axs[0].tick_params(axis="y",labelsize=11)
axs[0].tick_params(axis="x",labelsize=13)
axs[1].tick_params(axis="y",labelsize=11)
axs[1].tick_params(axis="x",labelsize=13)

#plt.ylabel(r"$\dot{E}/\epsilon$",fontsize=16)
#axs[1].legend(loc='lower right',title=r"$\kappa_{||} \sim 0$")
#plt.tight_layout()
plt.savefig("fcr_vs_pcpg_adiabatic_stackplot.pdf")
plt.close()



################################################
# Simulation details:
#   Lower Mach number: M ~ 0.15 instead of M ~ 0.5
#   CR transport now includes streaming
#   pc/pg = 1 = constant, now varying beta
#   "ad" denotes adiabatic eos instead of isothermal



beta = [1.5,7,50]
#fcr_iso_stream_res128 = [0.06, 0.16, 0.29]
fcr_ad_stream_res128 = [0.04, 0.13, 0.26]
fg_ad_stream_res128 = [0.87, 0.76, 0.64]
Heps_res128 = [0.41,0.67,0.55]

# res = 512, isothermal, with streaming
Heps_iso_res512_lowerMach = [0.596, 0.696]  # for beta = 1, 10
Heps_iso_res512_Mach05 = [0.33] # for beta = 1

Heps_iso_res256_lowerMach = [0.3788, 0.776, 0.723]

Heps_res256 = [0.40, 0.63, 0.54] # the beta = 1 outcome is much lower than the elongated adiabatic runs, why??

fcr_ad_stream_res128_elongated = [0.052,0.17,0.268]
fg_ad_stream_res128_elongated = [0.893,0.812,0.66]
Heps_res128_elongated = [0.71,0.67,0.55]

#vm = 10 instead of vm = 5
fcr_ad_stream_res128_elongated_vm10 = [0.06]
fg_ad_stream_res128_elongated_vm10 = [0.91]
Heps_res128_elongated_vm10 = [0.64]

fcr_ad_justStream_res128_elongated = [0.033, 0.094, 0.142]
fg_ad_justStream_res128_elongated = [0.908, 0.892,0.874]
Heps_justStream_res128_elongated = [0.7285,0.732,0.857]



plt.semilogx(beta,fcr_ad_stream_res128_elongated,'ko-',label=r"$f_{\rm CR}$",markersize=ms)
plt.semilogx(beta,fcr_ad_justStream_res128_elongated,'kd-',markersize=ms)
plt.semilogx(beta,np.array(fg_ad_stream_res128_elongated)-np.array(Heps_res128_elongated),'co-',label=r"$f_{\rm th} - f_{\rm CR,heating}$",markersize=ms)
plt.semilogx(beta,np.array(fg_ad_justStream_res128_elongated)-np.array(Heps_justStream_res128_elongated),'cd-',markersize=ms)
plt.semilogx(beta,fg_ad_stream_res128_elongated,'ro-',label=r"$f_{\rm th}$",markersize=ms)
plt.semilogx(beta,fg_ad_justStream_res128_elongated,'rd-',markersize=ms)
plt.ylim(0,1)
plt.xlim(1,60)
plt.xlabel(r"$\beta$",fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.text(1.5,0.5,r"$P_{CR}/P_{g} \sim 1$",fontsize=16,bbox=dict(edgecolor='black', alpha=0.1))
plt.title("With Streaming",fontsize=22)
#plt.legend(prop={'size': 12},ncol=3,bbox_to_anchor=(0.85, -0.2),frameon=False)
plt.legend(prop={'size': 12},ncol=2,bbox_to_anchor=(0.5, -0.12),title=r'$\bf{Color}$',   frameon=False)
add_legend([{'marker' : 'o', 'label' : r'$\kappa \sim 0.15 L_{0}v_{ph}$'}, {'marker' :   'd', 'label' : r'$\kappa = 0$'}],
          default_plot_args = {'ls' : '', 'c' : 'k', 'markersize' : '12','fillstyle' :    'full'}, ncol=1,bbox_to_anchor=(1.02, -0.12),frameon=False,
                                    title = r'$\bf{Symbol}$')
plt.tight_layout()
plt.savefig("fcr_vs_pcpg_stream.pdf")
plt.close()



# making a stacked and grouped bar chart plot with labels = no transport, diffusion only, streaming only, diffusion + streaming (isothermal), diffusion + streaming (adiabatic)

labels = ['$\beta \sim 1$', '$\beta \sim 10$', r'$\beta \sim 100$']
#arrays should be length 3, one for each beta

adiabatic_stream_diff_fheat = np.array([0.452, 0.812, 0.856])
adiabatic_stream_diff_fth = np.array([0.94, 0.89 ,0.89])
adiabatic_stream_diff_fcr = np.array([0.017, 0.054, 0.17])
adiabatic_streamOnly_fcr = np.array([0.016, 0.008, 0.029])
adiabatic_streamOnly_fth = np.array([0.9669, 1.007, 0.97])
adiabatic_streamOnly_fheat = np.array([0.456, 0.89, 0.93])

diff_fcr = [0.86, 0.87, 0.90]
diff_fE = [0.359, 0.27, 0.18]
diff_fheat = [0,0,0]
diff_fth = [0,0,0] # all zeros because no gas heating with isothermal eos

stream_fcr = [0.03, 0.036, 0.03]
stream_fE = [1.15, 0.44, 0.07]
stream_fheat = [0.386, 0.834, 0.89]
stream_fth = [0,0,0]

stream_diff_fcr = [0.03, 0.07, 0.15]
stream_diff_fE = [1.20, 0.396, 0.15]
stream_diff_fheat = [0.38, 0.827, 0.811]
stream_diff_fth = [0,0,0]

# I've created many more vectors than I planned to...let's put them in a DataFrame and
# write to a CSV file. In the future, I'll just add to that.

df_diff = pd.DataFrame([diff_fcr,diff_fheat,diff_fth],
        index = [r'f$_{\rm CR}$', r'f$_{\rm CR, heating}$', r'f$_{\rm th}$'],columns=[r'$\beta \sim 1$',r'$\beta \sim 10$', r'$\beta \sim   100$'])

df_stream = pd.DataFrame([stream_fcr,stream_fheat, stream_fth],
        index = [r'f$_{\rm CR}$', r'f$_{\rm CR, heating}$', r'f$_{\rm th}$'],columns=[r'$\beta \sim 1$',r'$\beta \sim 10$', r'$\beta \sim   100$'])

df_stream_diff = pd.DataFrame([stream_diff_fcr,stream_diff_fheat,stream_diff_fth],
        index = [r'f$_{\rm CR}$', r'f$_{\rm CR, heating}$', r'f$_{\rm th}$'],columns=[r'$\beta \sim 1$',r'$\beta \sim 10$', r'$\beta \sim   100$'])

df_adiabatic_stream_diff = pd.DataFrame([adiabatic_stream_diff_fcr,adiabatic_stream_diff_fheat,adiabatic_stream_diff_fth - adiabatic_stream_diff_fheat],
        index = [r'f$_{\rm CR}$', r'f$_{\rm CR, heating}$', r'f$_{\rm th}$'],columns=[r'$\beta \sim 1$',r'$\beta \sim 10$', r'$\beta \sim   100$'])

df_adiabatic_streamOnly = pd.DataFrame([adiabatic_streamOnly_fcr,adiabatic_streamOnly_fheat,adiabatic_streamOnly_fth - adiabatic_streamOnly_fheat],
        index = [r'f$_{\rm CR}$', r'f$_{\rm CR, heating}$', r'f$_{\rm th}$'],columns=[r'$\beta \sim 1$',r'$\beta \sim 10$', r'$\beta \sim   100$'])



# flatten into the right shape
df_flat = pd.concat([df_diff,df_stream,df_stream_diff,df_adiabatic_streamOnly, df_adiabatic_stream_diff])
print(df_flat.head(8))

# Write this stuff to a CSV file, so I can add to it easily
df_flat.to_csv('turbulent_energy_partitions.csv')


# this data is in wide form, but it's easier to work with in long form
df_flat = df_flat.to_numpy()
df_flat = df_flat.flatten('F')
print(df_flat)


# this works for now, since I haven't added more data to my CSV file. Need to change in future
dfall = pd.DataFrame(dict(
    beta = [r'$\beta \sim 1$',r'$\beta \sim 1$',r'$\beta \sim 1$',r'$\beta \sim 1$',r'$\beta \sim 1$',r'$\beta \sim 1$',r'$\beta \sim 1$',r'$\beta \sim 1$',r'$\beta \sim 1$',r'$\beta \sim 1$',r'$\beta \sim 1$',r'$\beta \sim 1$',r'$\beta \sim 1$',r'$\beta \sim 1$',r'$\beta \sim 1$',r'$\beta \sim 10$',
    r'$\beta \sim 10$',r'$\beta \sim 10$',r'$\beta \sim 10$',r'$\beta \sim 10$',r'$\beta \sim 10$',r'$\beta \sim 10$',r'$\beta \sim 10$',r'$\beta \sim 10$',r'$\beta \sim 10$',r'$\beta \sim 10$',r'$\beta \sim 10$',r'$\beta \sim 10$',r'$\beta \sim 10$',r'$\beta \sim 10$',
    r'$\beta \sim 100$',r'$\beta \sim 100$',r'$\beta \sim 100$',r'$\beta \sim 100$',r'$\beta \sim 100$',r'$\beta \sim 100$',r'$\beta \sim 100$',r'$\beta \sim 100$',r'$\beta \sim 100$',r'$\beta \sim 100$',r'$\beta \sim 100$',r'$\beta \sim 100$',r'$\beta \sim 100$',r'$\beta \sim 100$',r'$\beta \sim 100$'],
    labels = ["Diffusion Only","Diffusion Only","Diffusion Only","Streaming Only","Streaming Only","Streaming Only","Diff + Stream","Diff + Stream","Diff + Stream","Adiab. Stream Only","Adiab. Stream Only","Adiab. Stream Only","Adiab. Diff + Stream","Adiab. Diff + Stream","Adiab. Diff + Stream"] * 3,
    flabel = [r'$f_{\rm CR} \text{, CR Energization}$', r'$f_{\rm CR, heating} \text{, Streaming Energy Loss}$', r'$f_{\rm th} - f_{\rm CR, heating} \text{, Grid-Scale Heating}$'] * 15,
    f = df_flat))

print(dfall.head(16))

# Full
print("Full dataframe")
print(dfall)


##########################################################################3
# Now comes the fun part of plotting this as a stacked, grouped bar chart. After trying
# some methods in regular matplotlib, seaborn, etc., plotly emerged as the best option

# Without this ridiculous work-around, every PDF file I save has a box with "Loading MathJax" written in the lower left corner. Might be a problem with the Kaleido package used to write plotly images to files
figure="some_figure.pdf"
fig=px.scatter(x=[0, 1, 2, 3, 4], y=[0, 1, 4, 9, 16])
fig.write_image(figure, format="pdf")


time.sleep(2) # this delay gives some time for MathJax to load before the next figure (the real figure) is plotted and saved


fig = go.Figure()


# Note that all text that includes LaTeX anywhere near it has to be inside $ $
fig.update_layout(
    template="simple_white",
    xaxis=dict(ticklen=0),
    yaxis=dict(title_text=r"$\dot{E}/\epsilon$",range=(0,1.1)),
    font=dict(size=12),
    barmode="stack",
    legend=dict(
        font=dict(size=12),
        x=0.55,
        y=1.4,
        bgcolor='rgba(255, 255, 255, 0)',
        bordercolor='rgba(255, 255, 255, 0)'
    )
)

colors = ["blue", "red","orange"]

# Add traces for each bar. Grouped by $\beta \sim ...$, colored by f_{CR} or f_{CR, heating}
for r, c in zip(dfall.flabel.unique(), colors):
    plot_df = dfall[dfall.flabel == r]
    fig.add_trace(
        go.Bar(x=[plot_df.beta, plot_df.labels], y=plot_df.f, name=r, marker_color=c,opacity=0.7),
    )
fig.update_xaxes(
        tickangle = 90,
        tickson = "boundaries",
        ticks="inside",
        ticklen=0, # just get rid of the ticks
        dividerwidth=0,
        dividercolor='black',
        title_text = r"", # just get rid of the title
        title_font = {"size": 10},
        title_standoff = 20)

fig.write_image("barchart_plotly.pdf",engine="kaleido")


# Other variations on bar chart figures using Altair, seaborn. Keeping them here (commented out) in case they seem useful for the future.

"""

# Altair version -- mainly copied from stackexchange
def prep_df(df, name):
    df = df.stack().reset_index()
    df.columns = ['c1', 'c2', 'values']
    df['DF'] = name
    return df

df1 = prep_df(df_diff, 'DF1')
df2 = prep_df(df_stream, 'DF2')
df3 = prep_df(df_stream_diff, 'DF3')

df = pd.concat([df1, df2, df3])


alt.Chart(df).mark_bar().encode(

    # tell Altair which field to group columns on
    x=alt.X('c2:N', title=None),

    # tell Altair which field to use as Y values and how to calculate
    y=alt.Y('sum(values):Q',
        axis=alt.Axis(
            grid=False,
            title=None)),

    # tell Altair which field to use to use as the set of columns to be  represented in each group
    column=alt.Column('c1:N', title=None),

    # tell Altair which field to use for color segmentation
    color=alt.Color('DF:N',
            scale=alt.Scale(
                # make it look pretty with an enjoyable color pallet
                range=['#96ceb4', '#ffcc5c','#ff6f69'],
            ),
        ))\
    .configure_view(
        # remove grid lines around column clusters
        strokeOpacity=0
    )

"""

"""
# seaborn version -- mainly copied from stackexchange

df_diff['Name'] = "df1"
df_stream['Name'] = "df2"
df_stream_diff['Name'] = "df3"

dfall = pd.concat([pd.melt(i.reset_index(),
                           id_vars=["Name", "index"]) # transform in tidy format each df
                   for i in [df_diff, df_stream, df_stream_diff]],
                   ignore_index=True)

dfall.set_index(["Name", "index", "beta"], inplace=1)
dfall["vcs"] = dfall.groupby(level=["Name", "index"]).cumsum()
dfall.head(2)
dfall.reset_index(inplace=True)


c = ["blue", "green"]
for i, g in enumerate(dfall.groupby("beta")):
    ax = sns.barplot(data=g[1],
                     x="index",
                     y="vcs",
                     hue="Name",
                     color=c[i],
                     zorder=-i, # so first bars stay on top
                     edgecolor="k")
ax.savefig('barchart_test.pdf')
"""



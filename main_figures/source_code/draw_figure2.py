"""Render joint density, correlation strata and 18 complementary partitions."""
import csv
from pathlib import Path
import numpy as np
from scipy.stats import gaussian_kde, spearmanr
import matplotlib as mpl
mpl.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator

here = Path(__file__).resolve().parent
def read(name):
    with (here/"inputs/fig12"/f"{name}.csv").open() as stream:
        return list(csv.DictReader(stream))
navy,teal,orange,red="#1B4965","#2A9D8F","#D96B27","#B83A3A"
bg,ink,muted,grid="#F4F6F8","#1F2933","#52606D","#E2E8F0"
mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":7,"axes.labelsize":7,
    "xtick.labelsize":6.5,"ytick.labelsize":6.5,"legend.fontsize":6.5,
    "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none","axes.spines.top":False,
    "axes.spines.right":False,"axes.linewidth":.65,"xtick.major.width":.65,
    "ytick.major.width":.65,"xtick.major.size":2.2,"ytick.major.size":2.2,
    "axes.edgecolor":muted,"text.color":ink,"axes.labelcolor":ink,"xtick.color":muted,"ytick.color":muted})
W,H=167,186
fig=plt.figure(figsize=(W/25.4,H/25.4))
def axes(x,y,w,h):
    return fig.add_axes([x/W,y/H,w/W,h/H])
def text(x,y,s,**kw):
    fig.text(x/W,y/H,s,va="center",**kw)
for x,y,letter,title in [(3,181,"A","Treatment contrasts"),(85.5,181,"B","Between-study correlations"),
                         (3,89,"C","Within-study repeatability"),(85.5,89,"D","Treatment partitions")]:
    text(x,y,letter,fontsize=11,fontweight="bold",color="#1F2933")
    text(x+6.8,y,title,fontsize=8.5,fontweight="normal",color="#334155")

genes=read("Fig2A_all_9754_genes")
x,y=np.asarray([(float(r["root_delta"]),float(r["flower_delta"])) for r in genes]).T
rho=spearmanr(x,y).statistic
a=axes(15.7,111,55.3,46)
cmap=LinearSegmentedColormap.from_list("gene_density",["#E3EAF0","#81A8B9",navy])
hexes=a.hexbin(x,y,gridsize=39,mincnt=1,cmap=cmap,norm=LogNorm(),linewidths=.12,zorder=2)
a.margins(.025)
xlim=a.get_xlim()
ylim=a.get_ylim()
ref=dict(c=muted,lw=.5,ls=(0,(3,3)),zorder=2.5,alpha=.72)
a.axhline(0,**ref);a.axvline(0,**ref)
a.plot(xlim,xlim,**ref)
a.set(xlim=xlim,ylim=ylim,xlabel="PNS-H − Model, Δ mean log₂(1 + TPM)",
      ylabel="PNF − DM, Δ mean log₂(1 + TPM)")
a.text(.97,.96,f"n = {len(x):,}  |  ρ = {rho:.4f}",ha="right",va="top",
       transform=a.transAxes,fontsize=6.5,
       bbox=dict(boxstyle="round,pad=.28",fc="white",ec=grid,lw=.6,alpha=.93))
for values,limits,is_top in [(x,xlim,True),(y,ylim,False)]:
    margin=axes(15.7,157.6,55.3,7.0) if is_top else axes(71.7,111,6.5,46)
    coords=np.linspace(*limits,280)
    density=gaussian_kde(values)(coords)
    if is_top:
        margin.fill_between(coords,0,density,fc=navy,alpha=.19)
        margin.plot(coords,density,c=navy,lw=.8)
        margin.set(xlim=limits,ylim=(0,max(density)*1.06))
    else:
        margin.fill_betweenx(coords,0,density,fc=navy,alpha=.19)
        margin.plot(density,coords,c=navy,lw=.8)
        margin.set(ylim=limits,xlim=(0,max(density)*1.06))
    margin.axis("off")
text(16,172.3,"Genes per hexagon",fontsize=6.5,color=muted)
scale=axes(43,172.1,27,1.55)
cb=fig.colorbar(hexes,cax=scale,orientation="horizontal",ticks=[1,10,100])
cb.ax.set_xticklabels(["1","10","100"])
cb.solids.set_rasterized(False)
cb.ax.minorticks_off()
cb.ax.tick_params(length=1.5,pad=.8)

rows=read("cross_context_snr_calibration")
tpm="log2(1+TPM) mean difference"
items=[("treatment","TPM",tpm,"TPM"),("treatment","historical lengthScaledTPM","edgeR_logFC","logFC · A"),
       ("treatment","historical lengthScaledTPM","voom_t","voom t · A"),("disease","TPM",tpm,"TPM"),
       ("disease","root_A vs flower_A","edgeR_logFC","logFC · A"),("disease","root_A vs flower_A","voom_t","voom t · A"),
       ("disease","root_A vs flower_B2","edgeR_logFC","logFC · B2"),("disease","root_A vs flower_B2","voom_t","voom t · B2")]
b=axes(104.7,110.8,45.6,52.4)
bp=[8,7,6,4,3,2,1,0]
b.set(xlim=(-.112,.18),ylim=(-.65,8.85),yticks=bp,yticklabels=[i[3] for i in items],xlabel="Spearman ρ")
b.axhspan(5.25,8.85,fc=bg,zorder=-5)
b.axhspan(-.65,4.75,fc="#F3F8F7",zorder=-5)
b.axvline(0,c=muted,lw=.7,zorder=-1)
for loc in [-.1,.1]:b.axvline(loc,c=grid,lw=.4,zorder=-2)
b.set_xticks([-.1,0,.1])
b.tick_params(axis="y",length=0,pad=2)
b.spines["left"].set_visible(False)
text(153.4,166.7,"All ρ",fontsize=6.6,fontweight="bold")
for iy,item in zip(bp,items):
    vals={r["stratum"]:float(r["rho"]) for r in rows if
          (r["contrast"],r["count_scale"],r["statistic"])==item[:3]}
    b.plot([min(vals.values()),max(vals.values())],[iy,iy],c="#8699A6",lw=1,zorder=1)
    for stratum,col,mark in [("lower",orange,"D"),("higher",teal,"^"),("all",navy,"o")]:
        b.scatter(vals[stratum],iy,s=13,fc=col,marker=mark,ec="white",lw=.28,zorder=3)
    b.text(1.06,iy,f"{vals['all']:.4f}".replace("-","−"),transform=b.get_yaxis_transform(),
           va="center",color=navy,fontsize=6.5)
for yy,s in [(7,"Treatment"),(2,"Disease")]:
    b.text(-.4,yy,s,rotation=90,transform=b.get_yaxis_transform(),ha="center",va="center",fontsize=6.5,fontweight="bold")
handles=[Line2D([],[],marker=m,ls="none",mfc=cl,mec="none",markersize=4,label=name) for m,cl,name in
         [("D",orange,"Lower"),("o",navy,"All"),("^",teal,"Higher")]]
fig.legend(handles=handles,loc="center",bbox_to_anchor=(125/W,173/H),ncol=3,frameon=False,handlelength=1,columnspacing=1.1)

splits=read("flower_split_half_reliability")
c=axes(15.7,25.3,61.1,48.4)
all_rhos=[float(r["rho"]) for r in rows if r["stratum"]=="all" and r["rho"] not in ("NA","")]
low,high=round(min(all_rhos),3),round(max(all_rhos),3)
assert (low,high)==(-.035,.051)
c.axhspan(low,high,fc="#C6CED5",alpha=.58,zorder=-1)
c.set(ylim=(-.09,1.04),xlim=(.48,4.43),xticks=[1,2,3,4],xticklabels=["Treat.","Disease","Treat.","Disease"],
      ylabel="Split-half Spearman ρ",yticks=[0,.2,.4,.6,.8,1])
c.grid(axis="y",lw=.4,color=grid,zorder=-3)
for position,(stat,group,col) in enumerate([(tpm,"treatment",navy),(tpm,"disease",teal),
        ("voom_t","treatment",navy),("voom_t","disease",teal)],1):
    vals=np.asarray(sorted(float(r["rho"]) for r in splits if r["row_type"]=="pair" and r["contrast"]==group
                             and r["statistic"]==stat and r["scale"] in ("lengthScaledTPM","TPM")))
    assert len(vals)==18
    v=np.linspace(min(vals),max(vals),120)
    density=gaussian_kde(vals)(v)
    c.fill_betweenx(v,position-.04-density/max(density)*.38,position-.04,
                    fc=col,alpha=.32,ec=col,lw=.55,zorder=2)
    q1,q2,q3=np.percentile(vals,[25,50,75])
    c.bar(position,q3-q1,bottom=q1,width=.11,color="white",edgecolor=col,linewidth=.8,zorder=4)
    c.plot([position-.055,position+.055],[q2,q2],color=ink,lw=1,zorder=5)
    jitter=np.random.default_rng(141+position).uniform(.10,.29,len(vals))
    c.scatter(position+jitter,vals,s=6.8,fc=col,ec="none",alpha=.83,zorder=3)
    c.text(position,.977,f"{q2:.2f}",ha="center",va="center",fontsize=7.5,fontweight="bold",color=col)
for i,name in [(1.5,"Δ mean log TPM"),(3.5,"voom t · A")]:
    c.text(i,-.18,name,ha="center",transform=c.get_xaxis_transform(),fontsize=6.6)
# Key below the axes keeps the narrow reference band visible across the four groups.
fig.add_artist(Line2D([16/W,20.5/W],[9/H,9/H],transform=fig.transFigure,color="#BBC5CD",lw=3))
text(22,9,"Cross-study ρ span",fontsize=6.5,color=muted)

d=axes(97.7,25.3,65.5,48.4)
pairs=read("split_treatment_labels")
partners=["PNF19_4","PNF19_2","PNF19_1"]
colors=[navy,teal,orange]
d.set(xlim=(.45,18.55),ylim=(-.09,1.04),ylabel="Split-half Spearman ρ",xlabel="Balanced partition",
      xticks=[1,3,6,9,12,15,18],yticks=[0,.2,.4,.6,.8,1])
d.grid(axis="y",lw=.4,color=grid,zorder=-1)
for lane,(p,col) in enumerate(zip(partners,colors)):
    block=[r for r in pairs if r["PNF20_partner"]==p]
    ids=sorted({int(r["pair"]) for r in block})
    assert ids==list(range(1+6*lane,7+6*lane))
    d.axvspan(min(ids)-.5,max(ids)+.5,fc=col,alpha=.065,zorder=-3)
    d.text(np.mean(ids),1.035,p,ha="center",transform=d.get_xaxis_transform(),
           fontsize=6.65,fontweight="bold",color=col)
    for pair in ids:
        obs={r["statistic"]:float(r["rho"]) for r in block if int(r["pair"])==pair}
        d.plot([pair,pair],[obs[tpm],obs["voom_t"]],c=col,lw=.65)
        d.scatter(pair,obs[tpm],marker="o",fc=col,ec="white",lw=.2,s=13,zorder=3)
        d.scatter(pair,obs["voom_t"],marker="^",fc="white",ec=col,lw=.7,s=15,zorder=4)
    if lane==2:
        v=[float(r["rho"]) for r in block if r["statistic"]==tpm]
        span=f"{min(v):.3f} to {max(v):.3f}".replace("-","−")
        d.annotate("Pair 18",xy=(18,min(v)),xytext=(13.2,.115),ha="center",
            arrowprops=dict(arrowstyle="->",lw=.65,color=muted),
            bbox=dict(boxstyle="round,pad=.2",fc="white",ec=grid,lw=.55),fontsize=6.5)
        d.text(17.85,.95,"Fixed PNF, varying DM\nTPM: "+span,ha="right",va="top",
            bbox=dict(boxstyle="round,pad=.27",fc="white",ec=grid,lw=.55),fontsize=6.5,linespacing=1.4)
text(97.7,82,"PNF20_1 partner",fontsize=6.5,color=muted)
handles=[Line2D([],[],marker=m,ls="none",markersize=4,markeredgewidth=.8,mec=navy,mfc=fc,label=label)
         for m,fc,label in [("o",navy,"TPM"),("^","white","voom t")]]
fig.legend(handles=handles,loc="center",bbox_to_anchor=(126/W,9/H),ncol=2,frameon=False,handlelength=.9)

out=here.parent/"exports"
(out/"preview_png").mkdir(parents=True,exist_ok=True)
for ext in ("pdf","svg"):fig.savefig(out/f"Figure2.{ext}")
fig.savefig(out/"preview_png/Figure2.png",dpi=300)

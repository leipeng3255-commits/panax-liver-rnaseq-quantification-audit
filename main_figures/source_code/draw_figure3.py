"""Render Figure 3 from the archived gene and CAMERA tables. Run with Python 3."""
from pathlib import Path
import csv
import gzip
import math
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = Path(__file__).resolve().parent
INPUT = HERE / "inputs"
OUTPUT = HERE.parent / "exports"
W, H = 167, 215
NAVY, TEAL, ORANGE, RED = "#1B4965", "#2A9D8F", "#D96B27", "#B83A3A"
PALE, INK, GRAY, GRID = "#F4F6F8", "#1F2933", "#52606D", "#E2E8F0"
mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":7,"text.color":INK,
                     "axes.labelcolor":INK,"axes.edgecolor":GRAY,"xtick.color":GRAY,
                     "ytick.color":GRAY,"axes.linewidth":.65,"pdf.fonttype":42,
                     "ps.fonttype":42,"svg.fonttype":"none","axes.spines.top":False,
                     "axes.spines.right":False,"xtick.major.size":2.5,
                     "xtick.major.width":.65,"xtick.labelsize":6.5})
fig = plt.figure(figsize=(W/25.4, H/25.4))


def read(name):
    path = INPUT / name
    opener = gzip.open if name.endswith(".gz") else open
    with opener(path, "rt") as stream:
        return list(csv.DictReader(stream, delimiter="\t" if ".tsv" in name else ","))


def ax_at(x, y, w, h):
    return fig.add_axes([x/W, 1-(y+h)/H, w/W, h/H])


def text(x, y, label, **opts):
    fig.text(x/W, 1-y/H, label, va="center", **opts)


def title(x, y, letter, label):
    text(x, y, letter, fontsize=11,fontweight="bold",color="#1F2933")
    text(x+6.3, y, label, fontsize=8.5,fontweight="normal",color="#334155")


def key(x, y, marker, color, label, fill=None):
    fig.add_artist(plt.Line2D([x/W],[1-y/H],transform=fig.transFigure,marker=marker,
                             markersize=4,markeredgewidth=.75,markeredgecolor=color,
                             markerfacecolor=fill or color, linestyle="none"))
    text(x+2.4,y,label,fontsize=6.7)


def base_axis(ax):
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='x',pad=2)
    ax.set_axisbelow(True)


schemes = ("A","B1","B2","C")
colors = [NAVY,TEAL,ORANGE,RED]
shapes = ['o','s','^','D']
scale = {(r["scheme"],r["symbol"]):r for r in read("tximport_scale_sensitivity_13genes.csv")}
genes = sorted([g for s,g in scale if s=="A"], key=lambda g:float(scale["A",g]["edgeR_logFC"]),reverse=True)
libraries = read("ALL13_ALL12_GENE_SCALE_TRACE.tsv")
assert len(genes)==13 and genes[:2]==["Nat8l","Cdkn1a"]
minima={}
for g in genes:
    vals=[float(r["unscaled_estimated_fragments"]) for r in libraries if r["symbol"]==g and r["group"] in {"DM","PNF"}]
    assert len(vals)==8 and (sum(n>=10 for n in vals)>=4)==(scale["B2",g]["passed_filter"]=="TRUE")
    minima[g]=min(vals)
    assert math.isclose(minima[g],min(float(scale['A',g][j+'_unscaled_fragments_min']) for j in ('PNF','DM')),abs_tol=.001)
assert minima["Nat8l"]==0

title(4,5,"A","Flower treatment effects")
title(86,5,"B","Support and min. fragments")
for i,s in enumerate(schemes):
    key(9+i*15,12,shapes[i],colors[i],s)
effects=ax_at(22,21,60,72)
tiles=ax_at(88,21,37,72)
counts=ax_at(130,21,32,72)
for a in [effects,tiles,counts]:
    base_axis(a); a.set_ylim(12.6,-.6)
    for i in range(0,13,2):
        a.axhspan(i-.5,i+.5,color=PALE,linewidth=0,zorder=-2)
    a.axhline(1.5,color=GRAY,linestyle=(0,(2,2)),linewidth=.55,zorder=1)
effects.set_xlim(-3.15,6.25);effects.set_xticks([-2,0,2,4,6])
effects.grid(axis='x',color=GRID,linewidth=.4)
effects.axvline(0,color=GRAY,linewidth=.7)
effects.set_xlabel("edgeR log₂ fold change",labelpad=2.5)
tiles.set_xlim(-.5,3.5);tiles.axis("off")
for j,s in enumerate(schemes):
    text(88+(j+.5)*37/4,17.6,s,ha="center",fontweight="bold")
counts.set_xlim(-.09,2.87)
counts.set_xticks([math.log10(1+n) for n in (0,10,100,500)],["0","10","100","500"])
counts.axvline(math.log10(11),color=ORANGE,linestyle=(0,(3,2)),linewidth=.85,zorder=2)
text(146,12.2,"Min. unscaled fragments",ha="center",fontsize=6.7,fontweight="bold")
text(146,17.1,"log₁₀(1 + fragments)",ha="center",fontsize=6.4,color=GRAY)
states={"both":("D",NAVY,"white"),"ql":("o","#D7EBE8",TEAL),
        "ns":(".","#E8EDF0",GRAY),"excluded":("x","#F5E2DE",RED)}
for i,g in enumerate(genes):
    effects.text(-.035,i,g,ha="right",va="center",transform=effects.get_yaxis_transform(),
                 fontstyle="italic",fontsize=7.4,fontweight="bold" if i==0 else "normal")
    for j,s in enumerate(schemes):
        row=scale[s,g];y=i+[-.31,-.105,.105,.31][j]
        if row["passed_filter"]=="FALSE":
            state="excluded"
        else:
            fc=float(row["edgeR_logFC"])
            effects.plot([0,fc],[y,y],color=colors[j],linewidth=.45,alpha=.38)
            effects.plot(fc,y,marker=shapes[j],markersize=3.2,color=colors[j],markeredgewidth=.25,markeredgecolor="white")
            state="both" if float(row["QL_BH"])<.05 and float(row["TREAT_BH"])<.05 else ("ql" if float(row["QL_BH"])<.05 else "ns")
        mark,bg,fg=states[state]
        tiles.add_patch(FancyBboxPatch((j-.435,i-.41),.87,.82,boxstyle="round,pad=0.015,rounding_size=.14",facecolor=bg,edgecolor="white",linewidth=.55))
        tiles.plot(j,i,marker=mark,color=fg,markersize=3.7 if state!="ns" else 3,markeredgewidth=.85)
    n=minima[g]
    counts.barh(i,math.log10(1+n),height=.43,color="#7BA6AB",linewidth=0,zorder=1)
    if n==0:
        counts.plot(0,i,"o",color=RED,markersize=3)
        counts.annotate("0",(0,i),xytext=(4,0),textcoords="offset points",va="center",color=RED,fontsize=6.4)
delta=float(scale["B1","Nat8l"]["edgeR_logFC"])-float(scale["A","Nat8l"]["edgeR_logFC"])
text(24,17.8,"Nat8l",fontsize=6.6,fontstyle="italic",color=GRAY)
text(33,17.8,f"Δ = {delta:.2f} (A → B1)".replace("-","−"),fontsize=6.6,color=GRAY)
for x,st,label in [(89,"both","QL + TREAT"),(112,"ql","QL only"),(89,"ns","Tested NS"),(112,"excluded","Excluded")]:
    key(x,99 if st in ("both","ql") else 103.5,states[st][0],NAVY if st=="both" else states[st][2],label)
text(131,99,"Dashed line: 10 fragments",fontsize=6.2,color=GRAY)
text(131,103.5,"BH threshold: 0.05",fontsize=6.2,color=GRAY)

title(4,111,"C","CAMERA detections")
key(104,111,"o",NAVY,"cor = 0");key(128,111,"o",TEAL,"Estimated cor",fill="white")
overlap=read("camera_overlap.csv")
def treatment(tag):
    rows=read(tag+"_treatment_CAMERA.tsv.gz")
    rows=[r for r in rows if r["scenario"]=="full"]
    return [sum(float(r["BH"])<.05 for r in rows if r["correlation_setting"]==setting and r["BH"]!="NA") for setting in ("zero","estimated")]
def disease(tag,column):
    return [int(next(r for r in overlap if r['flower_scale']==tag and r['correlation']==k)[column]) for k in ("zero","estimated")]
contrasts=[("Root",treatment("root")),("Flower",treatment("flower")),("Root",disease("flower_A","root_hits")),
           ("Flower A",disease("flower_A","flower_hits")),("Flower B2",disease("flower_B2","flower_hits"))]
ca=ax_at(43,115.6,116,28);base_axis(ca);ca.set_ylim(4.75,-.75);ca.set_xlim(-.7,30)
ca.axhspan(1.5,4.65,color=PALE,linewidth=0,zorder=-2)
ca.set_xticks(range(0,31,5));ca.grid(axis="x",color=GRID,linewidth=.4)
for i,(label,(zero,estimated)) in enumerate(contrasts):
    ca.text(-.04,i,label,transform=ca.get_yaxis_transform(),ha="right",va="center",fontsize=7.1)
    ca.annotate("",(estimated+.28,i),(zero-.3,i),arrowprops={"arrowstyle":"->","lw":.85,"color":"#8BA7AD"})
    ca.scatter(zero,i,c=NAVY,s=16,zorder=3);ca.scatter(estimated,i,facecolors="white",edgecolors=TEAL,s=17,zorder=3,linewidth=.85)
    ca.annotate(str(zero),(zero,i),xytext=(5,0),textcoords="offset points",va="center",fontsize=7,color=NAVY,fontweight="bold")
    ca.annotate(str(estimated),(estimated,i),xytext=(-4.5,0),textcoords="offset points",ha="right",va="center",fontsize=6.8,color=TEAL)
for y,label in [(119,"Treatment"),(134.5,"Disease")]:
    text(4,y,label,color=GRAY,fontsize=7,fontweight="bold")
ca.set_xlabel("Significant Hallmark sets",labelpad=2)

title(4,157,"D","Shared disease Hallmarks (cor = 0)")
key(132,157,"o",NAVY,"Root");key(146,157,"D",ORANGE,"Flower A")
over=next(r for r in overlap if r["flower_scale"]=="flower_A" and r["correlation"]=="zero")
same=over["same_direction_names"].split("|");both=over["overlap_names"].split("|")
different=[k for k in both if k not in same]
assert len(same)==5 and len(different)==3
test={ (r["scale"],r["gene_set"]):r for r in read("disease_positive_control_camera.csv") if r["correlation"]=="zero"}
pathway_labels={"EPITHELIAL_MESENCHYMAL_TRANSITION":"Epithelial–mesenchymal transition","KRAS_SIGNALING_DN":"KRAS signaling DN","MYC_TARGETS_V1":"MYC targets V1"}
da=ax_at(56,161,105,44);base_axis(da);da.set_ylim(9.5,-.7);da.set_xlim(-20.5,20.5)
da.axhspan(5.6,9.45,color=PALE,linewidth=0,zorder=-2)
da.set_xticks([-20,-10,0,10,20]);da.grid(axis="x",color=GRID,linewidth=.4)
for q in [-1,1]:
    da.axvline(q*-math.log10(.05),c=GRAY,linewidth=.65,linestyle=(0,(2,2)))
da.axvline(0,c=GRAY,lw=.75)
da.text(2.3,0,"±1.30  ·  BH q = 0.05",color=GRAY,fontsize=6.2,va="center")
for group,y in [(same,0),(different,6)]:
    heading=('Concordant' if group is same else 'Discordant')+f' ({len(group)})'
    da.text(-.48,y,heading,transform=da.get_yaxis_transform(),fontsize=7.25,fontweight="bold",va='center')
    for idx,k in enumerate(group,1):
        yrow=y+idx
        vals=[]
        for dataset in ['root_A','flower_A']:
            r=test[dataset,k];assert float(r['BH'])<.05
            vals.append(-math.log10(float(r["BH"]))*(1 if r["direction"]=="Up" else -1))
        assert ((vals[0]*vals[1])>0)==(group is same)
        base=k.removeprefix("HALLMARK_");lab=pathway_labels.get(base,base.replace("_"," ").capitalize())
        da.text(-.016,yrow,lab,va="center",ha="right",transform=da.get_yaxis_transform(),fontsize=6.7)
        da.plot(vals,[yrow,yrow],c="#90A4AC",lw=.9)
        for v,c,m in zip(vals,[NAVY,ORANGE],['o','D']):
            da.plot(v,yrow,marker=m,ms=4,mec="white",mew=.35,mfc=c)
da.set_xlabel("Signed −log₁₀(BH q)   ·   downregulated  ←  0  →  upregulated",labelpad=2)
OUTPUT.mkdir(exist_ok=True)
(OUTPUT/"preview_png").mkdir(exist_ok=True)
for ext in ["pdf","svg"]:
    fig.savefig(OUTPUT/f"Figure3.{ext}")
fig.savefig(OUTPUT/"preview_png/Figure3.png",dpi=300)

"""Render study design and analysis topology at final publication size."""
from pathlib import Path
from collections import Counter
import csv
import gzip
import matplotlib as mpl
mpl.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, PathPatch, Circle, Ellipse
from matplotlib.path import Path as VectorPath
import numpy as np

here = Path(__file__).resolve().parent
inputs = here / "inputs"
def table(name):
    with (inputs/"fig12"/(name+".csv")).open() as stream:
        return list(csv.DictReader(stream))
libraries = Counter((r["study"],r["group"]) for r in table("library_groups"))
calibration=table("cross_context_snr_calibration")
design={r["metric"]:int(r["value"]) for r in table("design_counts")}
splits=table("split_treatment_labels")
genes=table("Fig2A_all_9754_genes")
with gzip.open(inputs/"flower_treatment_CAMERA.tsv.gz","rt") as stream:
    camera=[r for r in csv.DictReader(stream,delimiter="\t") if r["scenario"]=="full"]
n_sets=len({r["gene_set"] for r in camera})
domain={r["statistic"]:int(r["n_genes"]) for r in calibration if r["contrast"]=="treatment" and r["stratum"]=="all"}
mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":7,
                     "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})
navy,teal,orange,red="#1B4965","#2A9D8F","#D96B27","#B83A3A"
bg,ink,muted,grid="#F4F6F8","#1F2933","#52606D","#E2E8F0"
width,height=167,160
fig=plt.figure(figsize=(width/25.4,height/25.4))
ax=fig.add_axes([0,0,1,1],xlim=(0,width),ylim=(height,0),frameon=False)
ax.axis("off")
def txt(x,y,s,**kw):
    ax.text(x,y,s,va="center",color=kw.pop("color",ink),**kw)
def roundbox(x,y,w,h,color,edge="none",radius=2):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0,rounding_size={radius}",
                                 fc=color,ec=edge,lw=.65))
def stroke(x1,y1,x2,y2,col=muted,lw=.7):
    ax.plot([x1,x2],[y1,y2],color=col,lw=lw,solid_capstyle="round")
def arrow(x1,y1,x2,y2,col=muted):
    ax.annotate("",(x2,y2),(x1,y1),arrowprops=dict(arrowstyle="-|>",color=col,
                                                    lw=.8,shrinkA=0,shrinkB=0,mutation_scale=6))
def mouse(x,y,col):
    ax.add_patch(Ellipse((x+3.2,y),6,3.5,fc=col,ec="none"))
    ax.add_patch(Ellipse((x+6.2,y-.2),3.2,2.3,angle=10,fc=col,ec="none"))
    ax.add_patch(Circle((x+5.6,y-1.4),.85,fc=col))
    p=VectorPath([(x+.4,y+.3),(x-2.3,y+.2),(x-1.8,y-1.6),(x-3.0,y-1.3)],
                 [1,4,4,4])
    ax.add_patch(PathPatch(p,fc="none",ec=col,lw=.65))
def bracket(x,y1,y2,col):
    ax.plot([x,x+1.5,x+1.5],[y1,y1,y2],color=col,lw=1)
    arrow(x+1.5,y2,x-.5,y2,col)
def cohort(x,letter,study,title,color,mouse_desc,rows,contrasts):
    y,w,h=4,78.8,86.8
    roundbox(x,y,w,h,bg,grid)
    roundbox(x,y,w,11,color)
    ax.add_patch(Rectangle((x,y+6),w,5,fc=color,ec="none"))
    txt(x+2.6,y+5.6,letter,fontsize=11,fontweight="bold",color="white")
    txt(x+8.6,y+5.6,study+" · "+title,fontsize=9,fontweight="normal",color="white")
    mouse(x+5,y+16.5,color)
    txt(x+15,y+16.6,mouse_desc,fontsize=7.2)
    for j,(key,name,desc) in enumerate(rows):
        yy=y+26+j*11.1
        stroke(x+3,yy+5,x+68.5,yy+5,grid,.6)
        txt(x+3,yy,name,fontsize=7.7,fontweight="bold",color=color)
        txt(x+17,yy,desc,fontsize=6.5)
        txt(x+68.4,yy,f"n = {libraries[study,key]}",fontsize=6.5,ha="right")
    bracket(x+71,y+26,y+37.1,teal)
    bracket(x+74.6,y+37.1,y+48.2,orange)
    for i,(k,contrast,col) in enumerate([("Disease",contrasts[0],teal),("Treatment",contrasts[1],orange)]):
        yy=y+57+i*5.5
        stroke(x+3,yy,x+5.5,yy,col,1.7)
        txt(x+7,yy,k,fontsize=6.8,fontweight="bold")
        txt(x+28,yy,contrast,fontsize=6.8)
    txt(x+3,y+72.5,"Treatment week",fontsize=6.5,color=muted)
    x_left,x_right=x+30,x+73.8
    stroke(x_left,y+72.5,x_right,y+72.5,col=color)
    for tick in range(7):
        xx=x_left+(x_right-x_left)*tick/6
        stroke(xx,y+71.7,xx,y+73.7,col=color)
        if tick%2==0:txt(xx,y+77.2,str(tick),ha="center",fontsize=6.5)
cohort(3,"A","PRJNA1462623","Commercial saponins",navy,"Male ApoE-knockout mice",[
    ("Control","Control","Chow + saline"),("Model","Model","HFD + saline"),
    ("PNS","PNS-H","HFD + 225 mg/kg/d i.p. PNS-H")],["Model − Control","PNS-H − Model"])
cohort(85.2,"B","PRJNA1545393","Flower extract",teal,"Male BKS db/db mice and db/m littermates",[
    ("NC","NC","db/m + water"),("DM","DM","db/db + water"),
    ("PNF","PNF","db/db + 3.9 g/kg/d oral PNF")],["DM − NC","PNF − DM"])
txt(3,97,"C",fontsize=11,fontweight="bold",color="#1F2933")
txt(9,97,"Analytical workflow",fontsize=8.5,fontweight="normal",color="#334155")
stages=[("1. Within-study","contrasts",navy),("2. Cross-study","calibration",teal),
        ("3. Quantification","audit",orange),("4. Gene-set","dependence",navy)]
for i,(title,subtitle,color) in enumerate(stages):
    x=3+41.05*i
    roundbox(x,103,37.9,53,bg,grid)
    txt(x+2.2,109,title,fontweight="bold",fontsize=7.55,color=color)
    txt(x+2.2,113,subtitle,fontweight="bold",fontsize=7.55,color=color)
    if i<3:arrow(x+38.05,128.4,x+40.9,128.4,muted)
    # Compact vector pictograms use actual observations where axes imply data.
    if i==0:
        for j,n in enumerate([libraries["PRJNA1462623","PNS"],libraries["PRJNA1545393","PNF"]]):
            for k in range(n):
                ax.scatter(x+5+k*3,122+j*4,s=12,c=[navy,teal][j])
        arrow(x+19,124,x+25,124)
        txt(x+27,124,"Δ",fontsize=9,color=color)
        body=("GENCODE M37\nSalmon 2.8.0\nCoding genes\n"
              f"Root {design['root_treatment_genes']:,}\nFlower {design['flower_treatment_genes']:,}")
    elif i==1:
        pts=np.array([(float(g["root_delta"]),float(g["flower_delta"])) for g in genes[::270]])
        xx,yy=pts.T
        ax.scatter(x+9+2.1*xx,124-2.1*yy,s=2,fc=navy)
        stroke(x+3,129,x+15,129,grid)
        for j in range(4):ax.scatter(x+21+j*2.2,122,s=7,fc=teal if j<2 else navy)
        for j in range(4):ax.scatter(x+21+j*2.2,126,s=7,fc=navy if j<2 else teal)
        body=(f"{domain['log2(1+TPM) mean difference']:,} TPM genes\n{domain['edgeR_logFC']:,} model genes\n"
              f"{len(set(r['pair'] for r in splits))} balanced splits")
    elif i==2:
        for y,blocks in [(122,[(3,3),(10,3),(18,6)]),(126,[(14,3),(18,2)])]:
            stroke(x+blocks[0][0],y,x+sum(blocks[-1]),y,color)
            for start,length in blocks:ax.add_patch(Rectangle((x+start,y-1),length,2,fc=color))
        ax.add_patch(Circle((x+22.7,123.7),4,fc="none",ec=muted,lw=.8))
        stroke(x+25.6,126.6,x+28.1,129.1,muted,1.2)
        body=f"4 tximport schemes\n(A, B1, B2, C)\n{design['flower_transcripts']} transcripts"
    else:
        setting={s:sum(float(r["BH"])<.05 for r in camera if r["correlation_setting"]==s)
                 for s in set(r["correlation_setting"] for r in camera)}
        a,b=max(setting.values()),min(setting.values())
        stroke(x+4,128,x+29,128,grid)
        stroke(x+7,128-a*.55,x+23,128-b*.55,color,1)
        ax.scatter([x+7,x+23],[128-a*.55,128-b*.55],s=12,fc=color)
        body=f"{n_sets} Hallmark sets\nCAMERA cor = 0\nvs estimated"
    txt(x+2.2,142,body,fontsize=6.8,linespacing=1.45)
out=here.parent/"exports"
(out/"preview_png").mkdir(parents=True,exist_ok=True)
for ext in ("pdf","svg"):fig.savefig(out/f"Figure1.{ext}")
fig.savefig(out/"preview_png/Figure1.png",dpi=300)

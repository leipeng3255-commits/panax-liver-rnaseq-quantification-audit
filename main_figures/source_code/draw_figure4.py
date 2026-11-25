"""Render the Nat8l quantification audit from the archived model and transcript tables."""
import csv
from math import log2, log10
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Polygon
import numpy as np

HERE = Path(__file__).resolve().parent
def read(file, sep=","):
    with (HERE / "inputs" / file).open() as stream:
        return list(csv.DictReader(stream, delimiter=sep))

anno = {d["transcript"]: d for d in read("nat8l_annotation_structure.csv")}
gene = {d["isolate"]: d for d in read("ALL13_ALL12_GENE_SCALE_TRACE.tsv", "\t") if d["symbol"]=="Nat8l"}
tx = {(d["isolate"], d["transcript_name"]): d for d in read("NAT8L_ALL12_TWO_TRANSCRIPT_CONTRIBUTIONS.tsv","\t")}
models = {d["scheme"]: d for d in read("tximport_scale_sensitivity_13genes.csv") if d["symbol"]=="Nat8l"}
libs = sorted(k for k in gene if gene[k]["group"]=="NC") + sorted(k for k in gene if gene[k]["group"]=="DM") + ["PNF19_1","PNF19_2","PNF19_4","PNF20_1"]
assert len(libs)==len(set(libs))==len(gene)==12
navy,teal,orange,red="#1B4965","#2A9D8F","#D96B27","#B83A3A"
pale,ink,gray,grid="#F4F6F8","#1F2933","#52606D","#E2E8F0"
mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":7,"pdf.fonttype":42,"ps.fonttype":42,
                     "svg.fonttype":"none","svg.hashsalt":"Nat8l-figure4","hatch.linewidth":.4,
                     "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.65})
fig=plt.figure(figsize=(167/25.4,218/25.4),facecolor="white")
ax=fig.add_axes([0,0,1,1],xlim=(0,167),ylim=(218,0))
ax.set_axis_off()

def text(x,y,s,sz=7,color=ink,ha="left",**kw):
    return ax.text(x,y,s,va="center",ha=ha,fontsize=sz,color=color,zorder=8,**kw)
def rule(x,y,x2,y2,col=grid,lw=.4,**kw):
    return ax.plot([x,x2],[y,y2],color=col,lw=lw,solid_capstyle="butt",**kw)
def box(x,y,w,h,color,**kw):
    return ax.add_patch(Rectangle((x,y),w,h,facecolor=color,**kw))
def pill(x,y,w,h,bg):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0,rounding_size={h/2}",fc=bg,ec="none",zorder=1))
def header(y,panel,s):
    text(3.7,y,panel,11,"#1F2933",weight="bold");text(10,y,s,8.5,"#334155",weight="normal")
def arrow(x,y,x2,y2,col,lw=.7,scale=5):
    ax.add_patch(FancyArrowPatch((x,y),(x2,y2),arrowstyle="-|>",shrinkA=0,shrinkB=0,
                                 mutation_scale=scale,color=col,lw=lw,zorder=4))
def point(x,y,mark,col,size=17):
    ax.scatter(x,y,marker=mark,c=col,s=size,linewidths=.45,edgecolors="white",zorder=6)

# A. True genomic spans, followed by a linear zoom of the same two exon intervals.
header(4.6,"A","Transcript structure")
text(163,4.6,"GRCm39  ·  chr5 (+)",7,gray,"right")
long,short=anno["Nat8l-201"],anno["Nat8l-202"]
def exons(d): return [tuple(map(int,iv.split("-"))) for iv in d["exon_intervals"].split("|")]
lo,hi=int(long["genomic_start"]),int(long["genomic_end"])
def project(v):return 45+(v-lo)/(hi-lo)*115
ex1,ex2=exons(short)
e1=ex1[1]-ex1[0]+1;e2=ex2[1]-ex2[0]+1
assert e1+e2==int(short["full_length_nt"]) and e2==int(short["exonic_overlap_with_201_nt"])
zoom_lo=(ex1[0]//1000)*1000;zoom_hi=((ex2[1]+499)//500)*500
x0,x1=project(zoom_lo),project(zoom_hi)
box(x0,9,x1-x0,15.7,"#FFF2E8",ec=orange,lw=.65,ls=(0,(3,2)),zorder=1)
for datum,base,color in [(long,13.2,navy),(short,21.2,orange)]:
    spans=exons(datum)
    rule(project(spans[0][0]),base,project(spans[-1][1]),base,color,.7,zorder=3)
    for start,stop in spans:
        box(project(start),base-1.25,project(stop+1)-project(start),2.5,color,ec="none",zorder=4)
    for (_,last),(nxt,_) in zip(spans[:-1],spans[1:]):
        for mid in np.linspace(project(last),project(nxt),4)[1:-1]:
            rule(mid-.45,base-.6,mid+.3,base,color,.55)
            rule(mid+.3,base,mid-.45,base+.6,color,.55)
    text(4,base,datum["transcript"],7.7,color,style="italic",weight="bold")
    if datum is short:
        text(4,base+3.8,f'{int(datum["full_length_nt"])} nt · CDS not defined',6.8,gray)
        continue
    text(23,base,f'{int(datum["full_length_nt"]):,} nt',6.8,gray)
text(45,9.4,"3 exons · protein coding",6.5,gray)
bar=project(lo+1000)-project(lo);rule(160-bar,20.5,160,20.5,ink,.8);text(160-bar/2,24,"1 kb",6.5,gray,"center")
ax.add_patch(Polygon([(x0,24.8),(x1,24.8),(159,33),(29,33)],closed=True,fc="#F9E4D7",
                     ec="none",alpha=.62,zorder=.7))
rule(x0,24.8,29,33,orange,.45,alpha=.45);rule(x1,24.8,159,33,orange,.45,alpha=.45)
def zoom(v):return 29+(v-zoom_lo)/(zoom_hi-zoom_lo)*130
text(4,35.9,"Zoom",7.2,gray,weight="bold")
text(29,34.8,f'{zoom_lo/1000:,.0f} kb',6.5,gray)
text(159,34.8,f'{zoom_hi/1000:,.1f} kb',6.5,gray,"right")
ref_y,zbase=37.7,42.1
exon3_start,exon3_end=exons(long)[2]
shared_start,shared_end=max(exon3_start,ex2[0]),min(exon3_end,ex2[1])
assert (shared_start,shared_end)==(34158124,34158187) and shared_end-shared_start+1==64
box(zoom(shared_start),ref_y-1.8,zoom(shared_end+1)-zoom(shared_start),zbase-ref_y+3.6,"#AAB4BE",alpha=.38,ec="none",zorder=1)
rule(zoom(zoom_lo),ref_y,zoom(exon3_start),ref_y,navy,.7,zorder=3)
for mid in (zoom_lo+320,zoom_lo+700):
    rule(zoom(mid)-.6,ref_y-.65,zoom(mid)+.35,ref_y,navy,.6,zorder=4)
    rule(zoom(mid)+.35,ref_y,zoom(mid)-.6,ref_y+.65,navy,.6,zorder=4)
box(zoom(exon3_start),ref_y-1,zoom(zoom_hi)-zoom(exon3_start),2,navy,ec="none",zorder=4)
text(27,ref_y,"201",6.6,navy,"right")
text(27,zbase,"202",6.6,orange,"right")
rule(zoom(ex1[1]),zbase,zoom(ex2[0]),zbase,orange,.85)
for start,stop in (ex1,ex2):
    box(zoom(start),zbase-1.65,zoom(stop+1)-zoom(start),3.3,orange,ec="none",zorder=4)
text((zoom(ex1[1])+zoom(ex2[0]))/2,zbase-1.9,"Intron",6.6,gray,"center")
text(zoom(ex1[0])-2,45.4,f'Exon 1 · {e1} nt',7,ink,weight="bold")
text(zoom(ex1[0])-2,49.1,"Within a 201 intron",6.7,gray)
text(zoom(ex2[0])-5,45.4,f'Exon 2 · {e2} nt',7,ink,weight="bold")
text(zoom(ex2[0])-5,49.1,"Shared with 201 exon 3",6.7,gray)
eff=[float(tx[(k,"Nat8l-202")]["EffectiveLength"]) for k in libs]

pill(56,51.9,54,5.3,"#FBEADE")
text(83,54.5,f"Effective length ≈ {min(eff):.0f}–{max(eff):.0f} nt",6.8,orange,"center",weight="bold")

# B. Aligned per-library tracks; a zero count does not define a TPM share.
rule(3.7,58.5,163,58.5)
header(63,"B","Library transcript shares")
for x,t,c in [(128,"201",navy),(148,"202",orange)]:
    box(x,61.8,2.4,2.4,c,ec="none");text(x+3.6,63,t,7)
xp=np.linspace(28,158,12)
for j,g in enumerate(["NC","DM","PNF"]):
    xx=xp[4*j:4*j+4].mean();pill(xx-18.8,67.3,37.6,4.6,pale);text(xx,69.6,g,7,gray,"center",weight="bold")
def skeleton(zero,height,maximum,ticks,pct=False):
    for val in ticks:
        yy=zero-height*val/maximum
        rule(21,yy,163,yy,grid if val else gray,.4 if val else .65,zorder=.8)
        text(19.8,yy,f"{val}%" if pct else str(val),6.5,gray,"right")
skeleton(90,13.5,12,[0,4,8,12]);skeleton(110.5,14.2,100,[0,50,100],True)
text(4.6,81.4,"Estimated\nfragments",6.6,gray,"center",rotation=90)
text(4.6,104.2,"Gene TPM\nshare",6.6,gray,"center",rotation=90)
pnf_shares=[]
for i,(x,k) in enumerate(zip(xp,libs)):
    rows=[tx[k,f"Nat8l-{n}"] for n in [201,202]]
    y=90
    for r,c in zip(rows,[navy,orange]):
        h=float(r["NumReads_estimated_fragments"])*13.5/12
        box(x-2.35,y-h,4.7,h,c,ec="none");y-=h
    if k.startswith("PNF"):
        text(x,y-1.9,f'{float(rows[1]["NumReads_estimated_fragments"]):.1f}',6.5,orange,"center",weight="bold")
    if float(rows[0]["gene_TPM"])>0:
        y=110.5
        for r,c in zip(rows,[navy,orange]):
            h=float(r["fraction_gene_TPM"])*14.2;box(x-2.35,y-h,4.7,h,c,ec="none");y-=h
        if k.startswith("PNF"):pnf_shares.append(float(rows[1]["fraction_gene_TPM"])*100)
    else:
        text(x,104,"n/a",6.5,gray,"center")
    text(x,113.9,k[len(gene[k]["group"]):],6.5,gray,"center")
span=f"{min(pnf_shares):.2f}–{max(pnf_shares):.2f}%"
text(xp[-4:].mean(),93.3,span+" from 202",6.5,orange,"center",weight="bold")
text(4,119,f"PNF: approximately one estimated 202 fragment per library contributed {span} of gene TPM.",6.8)
text(4,123.4,"n/a: zero total gene TPM. Library identifiers pair the group prefix with the suffix shown.",6.5,gray)

# C. Both measures refer to the same library. Lines have small horizontal offsets for legibility.
rule(3.7,127.4,163,127.4)
header(132,"C","Count amplification")
treat=libs[4:];xps=np.linspace(20,154,8)
raw=np.array([float(gene[k]["unscaled_estimated_fragments"]) for k in treat])
scaled=np.array([float(gene[k]["reconstructed_lengthScaledTPM_count"]) for k in treat])
badge=f"{(raw>=10).sum()} of 8 libraries ≥10 unscaled vs {(scaled>=10).sum()} of 8 length-scaled"
pill(86,129.4,77,5.3,pale);text(124.5,132,badge,6.5,gray,"center",weight="bold")
yc=lambda v:168.4-v*.45
box(13,yc(10),150,yc(0)-yc(10),"#F9F0E9",ec="none")
for v in [0,10,30,50]:
    rule(13,yc(v),163,yc(v),gray if v==10 else grid,.65 if v==10 else .4,ls=(0,(3,2)) if v==10 else "-")
    text(11.6,yc(v),v,6.5,gray,"right")
text(4.7,157,"Counts",6.6,gray,"center",rotation=90)
for i,(k,x,r,s) in enumerate(zip(treat,xps,raw,scaled)):
    if i>=4:
        arrow(x-1.5,yc(r),x+1.38,yc(s)+1.25,teal,.7,6)
        text(x+2,yc(s)-2.7,f"{r:.0f}→{s:.0f}",6.6,teal,"center",weight="bold")
    else:rule(x-1.5,yc(r),x+1.5,yc(s),gray,.75)
    point(x-1.5,yc(r),"o",navy);point(x+1.5,yc(s),"^",teal,20)
    text(x,171.4,k,6.5,gray,"center")
point(16,138.7,"o",navy);text(18.3,138.7,"Unscaled fragments",6.6)
point(47,138.7,"^",teal);text(49.3,138.7,"lengthScaledTPM",6.6)
text(17,161,"Below count ≥ 10 filter",6.6,gray,weight="normal")
rule(45,161.3,50,164.9,gray,.5)

# D. Numeric axes retain measured coefficients and BH-adjusted significance.
rule(3.7,175.7,163,175.7)
text(3.7,180.5,"D",11,"#1F2933",weight="bold")
text(10,180.5,"Nat8l",8.5,"#334155",style="italic",weight="normal")
text(20.5,180.5,"threshold support",8.5,"#334155",weight="normal")
starts=[35,106];widths=[56,54]
texts=["edgeR log₂ fold change","−log₁₀(TREAT q)"]
thresholds=[log2(1.5),-log10(.05)]
for j,(x,w,title,thr) in enumerate(zip(starts,widths,texts,thresholds)):
    text(x,186,title,7.3,weight="bold")
    refx=x+w*thr/6.6
    rule(refx,188,refx,204.5,gray,.6,ls=(0,(2.5,2)))
    for v in (0,2,4,6):
        xval=x+w*v/6.6;rule(xval,188,xval,204.4,grid,.4)
        text(xval,207,v,6.5,gray,"center")
    rule(x,204.4,x+w,204.4,gray,.65)
    text(x,215.3,"Reference: log₂(1.5) = 0.585" if j==0 else "Reference: −log₁₀(0.05) = 1.30",6.5,gray)
for i,(s,name) in enumerate([("A","lengthScaledTPM"),("C","scaledTPM"),("B1","fixed-domain offset")]):
    r=models[s];yy=190+i*5.5
    text(4,yy,s,7.1,weight="bold");text(9,yy,name,6.5,gray)
    for x,w,val,c in zip(starts,widths,[float(r["edgeR_logFC"]),-log10(float(r["TREAT_BH"]))],[navy,teal]):
        xv=x+w*val/6.6;rule(x,yy,xv,yy,c,1.5);point(xv,yy,"o",c,17);text(xv+1.9,yy,f"{val:.2f}",6.6,c)
assert models["B2"]["passed_filter"]=="FALSE"
text(4,211.1,"B2",7.1,weight="bold");text(9,211.1,"raw-filter offset",6.5,gray)
box(35,209.2,125,3.8,pale,edgecolor="#B8C0C7",hatch="////",linewidth=.3)
text(97.5,211.1,"Excluded before testing (raw filter)",6.6,gray,"center",
     bbox={"facecolor":pale,"edgecolor":"none","pad":1.2})
out=HERE.parent/"exports";(out/"preview_png").mkdir(parents=True,exist_ok=True)
for ext in ["pdf","svg","png"]:
    where=out/"preview_png" if ext=="png" else out
    fig.savefig(where/f"Figure4.{ext}",dpi=300)

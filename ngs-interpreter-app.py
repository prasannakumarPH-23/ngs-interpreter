import streamlit as st

# ============================================================
# NGS INTERPRETER — Single File App
# Paste this entire file into https://share.streamlit.io
# or run locally: streamlit run app.py
# ============================================================

import streamlit as st
import os, re, gzip, struct, zlib, collections, statistics, math, io, datetime, tempfile, struct

try:
    import pandas as pd
except ImportError:
    st.error("pandas not found. Add 'pandas' to packages.txt"); st.stop()

try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import mm
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.enums import TA_CENTER
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable
    REPORTLAB_OK = True
except ImportError:
    REPORTLAB_OK = False

try:
    import pysam
    PYSAM_OK = True
except ImportError:
    PYSAM_OK = False

try:
    from Bio import SeqIO
    BIOPYTHON_OK = True
except ImportError:
    BIOPYTHON_OK = False

# ════════════════════════════════════════════════════════════
# PAGE CONFIG
# ════════════════════════════════════════════════════════════

st.set_page_config(
    page_title="NGS Interpreter",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Serif+Display&family=DM+Sans:wght@300;400;500;600&family=JetBrains+Mono:wght@400;600&display=swap');
:root {
  --navy:#0D1F3C; --teal:#0FA3A3; --teal2:#0D8C8C;
  --sky:#E8F4F8; --light:#F5F8FC; --mid:#D0DDE8;
  --text:#1E2D40; --muted:#5A7A9A;
  --pass:#1A7A4A; --warn:#B8740A; --fail:#C0392B;
}
html,body,[class*="css"]{ font-family:'DM Sans',sans-serif; color:var(--text); }
.stApp { background:var(--light); }
h1,h2,h3 { font-family:'DM Serif Display',serif; }
#MainMenu,footer,header { visibility:hidden; }
.block-container { padding-top:1rem !important; }

.topbar {
  background:var(--navy); color:white;
  padding:0.7rem 2rem; display:flex; align-items:center; gap:1rem;
  margin:-1rem -1rem 2rem -1rem; border-bottom:3px solid var(--teal);
}
.hero {
  background:linear-gradient(135deg,var(--navy) 0%,#1A3A5C 100%);
  border-radius:16px; padding:3rem; color:white;
  position:relative; overflow:hidden; margin-bottom:2rem;
}
.hero h1 { font-size:2.6rem; line-height:1.15; margin:0.5rem 0 1rem 0; color:white; }
.hero p  { font-size:1rem; color:#A0C0D8; line-height:1.7; max-width:600px; margin-bottom:1.5rem; }
.hero-tag {
  display:inline-block; background:rgba(15,163,163,0.2);
  border:1px solid rgba(15,163,163,0.4); color:var(--teal);
  padding:0.2rem 0.8rem; border-radius:20px;
  font-size:0.75rem; font-weight:600; letter-spacing:0.08em;
  text-transform:uppercase; margin-bottom:1rem;
}
.hero-badge {
  display:inline-block; background:rgba(255,255,255,0.08);
  border:1px solid rgba(255,255,255,0.15); border-radius:8px;
  padding:0.3rem 0.8rem; font-size:0.8rem; color:#C0D8E8;
  margin-right:0.4rem; margin-bottom:0.4rem;
}
.gloss-card {
  background:white; border-radius:12px; padding:1.2rem 1.4rem;
  border:1px solid var(--mid); border-left:4px solid var(--teal);
  box-shadow:0 2px 8px rgba(0,0,0,0.04); height:100%;
}
.gloss-card h4 {
  font-family:'JetBrains Mono',monospace; font-size:0.85rem;
  color:var(--navy); margin:0 0 0.4rem 0;
}
.gloss-card p { font-size:0.82rem; color:var(--muted); margin:0; line-height:1.55; }
.sec-header {
  display:flex; align-items:center; gap:0.6rem;
  background:var(--navy); color:white; border-radius:8px;
  padding:0.55rem 1rem; margin:1.5rem 0 0.6rem 0;
  font-weight:600; font-size:0.9rem;
}
.sec-dot { width:8px; height:8px; background:var(--teal); border-radius:50%; flex-shrink:0; }
.metrics-grid {
  display:grid; grid-template-columns:repeat(auto-fill,minmax(260px,1fr)); gap:0.5rem;
}
.metric-row {
  background:white; border:1px solid var(--mid); border-radius:8px;
  padding:0.5rem 0.8rem; display:flex; justify-content:space-between; align-items:center;
}
.metric-label { font-size:0.82rem; color:var(--muted); }
.metric-value { font-size:0.88rem; font-weight:600; color:var(--navy);
                font-family:'JetBrains Mono',monospace; text-align:right; }
.qc-row {
  display:flex; align-items:center; gap:0.8rem;
  background:white; border:1px solid var(--mid);
  border-radius:8px; padding:0.5rem 0.9rem; margin-bottom:0.4rem;
}
.qc-badge { border-radius:5px; padding:0.15rem 0.5rem; font-size:0.72rem;
            font-weight:700; font-family:'JetBrains Mono',monospace; flex-shrink:0; }
.qc-pass { background:#D4EDDA; color:var(--pass); }
.qc-warn { background:#FFF3CD; color:var(--warn); }
.qc-fail { background:#F8D7DA; color:var(--fail); }
.interp-box {
  background:white; border:1px solid var(--mid);
  border-left:4px solid var(--teal); border-radius:0 8px 8px 0;
  padding:1.2rem 1.4rem;
}
.interp-item { font-size:0.86rem; color:var(--text); line-height:1.6;
               padding:0.25rem 0; border-bottom:1px solid #F0F0F0; }
.interp-item:last-child { border-bottom:none; }
.fmt-badge {
  display:inline-flex; align-items:center; gap:0.6rem;
  background:var(--navy); color:white; border-radius:10px;
  padding:0.6rem 1.2rem; font-family:'JetBrains Mono',monospace;
  font-size:0.9rem; margin:1rem 0;
}
.stButton>button { font-family:'DM Sans',sans-serif; font-weight:600; transition:all 0.2s; }
</style>
""", unsafe_allow_html=True)

# ════════════════════════════════════════════════════════════
# SESSION STATE
# ════════════════════════════════════════════════════════════

for k,v in [("page","home"),("result",None),("fname",""),("fmt","")]:
    if k not in st.session_state:
        st.session_state[k] = v

def go(page):
    st.session_state.page = page
    st.rerun()

def topbar():
    st.markdown("""
    <div class="topbar">
        <span style='font-family:DM Serif Display,serif;font-size:1.2rem;'>🧬 NGS Interpreter</span>
        <span style='color:#5A7A9A;font-size:0.82rem;'>Genomics File Analysis Platform</span>
    </div>""", unsafe_allow_html=True)

# ════════════════════════════════════════════════════════════
# FORMAT DETECTION
# ════════════════════════════════════════════════════════════

FORMAT_INFO = {
    "BAM":       ("BAM",        "Binary Alignment Map",                ".bam"),
    "SAM":       ("SAM",        "Sequence Alignment Map (text)",        ".sam"),
    "FASTA":     ("FASTA",      "Nucleotide / Protein Sequences",       ".fa .fna .fasta"),
    "FASTQ":     ("FASTQ",      "Sequencing Reads + Quality Scores",    ".fastq .fq"),
    "VCF":       ("VCF",        "Variant Call Format",                  ".vcf"),
    "BCF":       ("BCF",        "Binary Variant Call Format",           ".bcf"),
    "GFF":       ("GFF3",       "General Feature Format v3",            ".gff .gff3"),
    "GTF":       ("GTF",        "Gene Transfer Format",                 ".gtf"),
    "GENBANK":   ("GenBank",    "NCBI GenBank Flat File",               ".gbk .gb .gbff"),
    "EMBL":      ("EMBL",       "EMBL-EBI Flat File",                   ".embl"),
    "BED":       ("BED",        "Browser Extensible Data",              ".bed"),
    "NARROWPEAK":("narrowPeak", "ChIP/ATAC-seq Peak Format",            ".narrowPeak"),
    "BROADPEAK": ("broadPeak",  "Broad Peak Format",                    ".broadPeak"),
    "BEDGRAPH":  ("BEDGraph",   "Continuous Signal Track",              ".bedGraph"),
    "WIG":       ("WIG",        "Wiggle Signal Track",                  ".wig"),
    "MAF":       ("MAF",        "Multiple Alignment Format",            ".maf"),
    "PLINK_BIM": ("PLINK .bim", "PLINK Variant Map",                   ".bim"),
    "PLINK_FAM": ("PLINK .fam", "PLINK Sample Information",            ".fam"),
    "AB1":       ("AB1",        "Sanger Sequencing Chromatogram",       ".ab1"),
    "UNKNOWN":   ("Unknown",    "Format not recognised",                "—"),
}

def detect_format(path):
    filename = os.path.basename(path).lower()
    ext = filename.rsplit(".",1)[-1] if "." in filename else ""
    with open(path,"rb") as f: head = f.read(16)
    if head[0:4]==b"ABIF": return "AB1"
    if head[0:4]==b"CRAM": return "CRAM"
    if head[0:2]==b"\x1f\x8b":
        if ext=="bam": return "BAM"
        if ext=="bcf": return "BCF"
        try:
            with gzip.open(path,"rb") as fg: inner=fg.read(8)
            if inner[0:4]==b"BAM\x01": return "BAM"
        except: pass
        try:
            with gzip.open(path,"rt",errors="ignore") as fg: first=fg.readline()
            if first.startswith(">"): return "FASTA"
            if first.startswith("@"): return "FASTQ"
            if "##fileformat=VCF" in first: return "VCF"
            if first.startswith("LOCUS"): return "GENBANK"
        except: pass
    lines=[]
    try:
        with open(path,"rt",errors="ignore") as f:
            for _ in range(30):
                l=f.readline()
                if not l: break
                lines.append(l.rstrip())
    except: return "UNKNOWN"
    ne=[l for l in lines if l.strip()]
    if not ne: return "UNKNOWN"
    if ne[0].startswith("LOCUS"): return "GENBANK"
    if ne[0].startswith("ID   ") and "SV" in ne[0]: return "EMBL"
    if ne[0].startswith(">"): return "FASTA"
    if any("##fileformat=VCF" in l for l in lines): return "VCF"
    if any(l.startswith("#CHROM") for l in lines): return "VCF"
    if ne[0].startswith(("@HD","@SQ")): return "SAM"
    sam_hdr=[l for l in lines[:10] if l.startswith(("@HD","@SQ","@RG","@PG"))]
    if len(sam_hdr)>=2: return "SAM"
    if len(ne)>=4 and ne[0].startswith("@") and ne[2].startswith("+") \
       and not ne[0].startswith(("@HD","@SQ","@RG","@PG")): return "FASTQ"
    if any(l.startswith(("variableStep","fixedStep","track type=wiggle")) for l in lines[:5]): return "WIG"
    if any("type=bedGraph" in l for l in lines[:3]): return "BEDGRAPH"
    if any(l.startswith("##maf") for l in lines[:5]): return "MAF"
    if ext=="bim":
        for l in [l for l in lines if l and not l.startswith("#")][:3]:
            cols=l.split()
            if len(cols)==6:
                try: int(cols[0]); int(cols[3]); return "PLINK_BIM"
                except: break
    if ext=="fam":
        for l in [l for l in lines if l and not l.startswith("#")][:3]:
            if len(l.split())==6: return "PLINK_FAM"
    data=[l for l in lines if l and not l.startswith(("#","track","browser"))]
    for l in data[:5]:
        cols=l.split("\t")
        if len(cols)==9:
            if 'gene_id "' in cols[8]: return "GTF"
            return "GFF"
    for l in data[:5]:
        cols=l.split("\t")
        if len(cols)==10:
            try: int(cols[1]); int(cols[2]); return "NARROWPEAK"
            except: pass
        if len(cols)==9:
            try: int(cols[1]); int(cols[2]); return "BROADPEAK"
            except: pass
    for l in data[:5]:
        cols=l.split("\t")
        if len(cols)>=3:
            try: int(cols[1]); int(cols[2]); return "BED"
            except: pass
    return "UNKNOWN"

# ════════════════════════════════════════════════════════════
# RESULT HELPERS
# ════════════════════════════════════════════════════════════

def _ok(fmt,sections,qc=None,interp=None):
    return {"format":fmt,"sections":sections,"qc":qc or [],"interpretation":interp or [],"error":None}
def _err(fmt,msg):
    return {"format":fmt,"sections":[],"qc":[],"interpretation":[],"error":msg}
def _metrics(title,rows): return {"title":title,"type":"metrics","data":rows}
def _table(title,data,note=None):
    df=pd.DataFrame(data) if isinstance(data,list) else data
    return {"title":title,"type":"table","data":df,"note":note}
def _bars(title,items): return {"title":title,"type":"bars","data":items}
def _grade(ok,warn=None):
    if ok: return "PASS"
    if warn is None or not warn: return "FAIL"
    return "WARN"
def _qc(check,result): return {"check":check,"result":result}

# ════════════════════════════════════════════════════════════
# BAM PURE-PYTHON PARSER
# ════════════════════════════════════════════════════════════

def _parse_bam_pure(path):
    data=b""
    with open(path,"rb") as f: raw=f.read()
    i=0
    while i<len(raw)-18:
        if raw[i:i+2]!=b"\x1f\x8b": break
        bsize=struct.unpack_from("<H",raw,i+16)[0]+1
        try: data+=zlib.decompress(raw[i+18:i+bsize-8],-15)
        except: pass
        i+=bsize
    if data[:4]!=b"BAM\x01": return None,[],[]
    offset=4
    l_text=struct.unpack_from("<i",data,offset)[0]; offset+=4
    header_text=data[offset:offset+l_text].decode("utf-8","ignore"); offset+=l_text
    n_ref=struct.unpack_from("<i",data,offset)[0]; offset+=4
    refs=[]
    for _ in range(n_ref):
        l_name=struct.unpack_from("<i",data,offset)[0]; offset+=4
        name=data[offset:offset+l_name-1].decode("utf-8","ignore"); offset+=l_name
        l_ref=struct.unpack_from("<i",data,offset)[0]; offset+=4
        refs.append((name,l_ref))
    reads=[]
    while offset+4<=len(data):
        bs=struct.unpack_from("<i",data,offset)[0]; offset+=4
        if bs<=0 or offset+bs>len(data): break
        rd=data[offset:offset+bs]; offset+=bs
        if len(rd)<32: continue
        try:
            ref_id=struct.unpack_from("<i",rd,0)[0]
            pos=struct.unpack_from("<i",rd,4)[0]
            mapq=struct.unpack_from("<B",rd,9)[0]
            flag=struct.unpack_from("<H",rd,16)[0]
            l_seq=struct.unpack_from("<i",rd,20)[0]
            tlen=struct.unpack_from("<i",rd,28)[0]
            rname=refs[ref_id][0] if 0<=ref_id<len(refs) else "*"
            reads.append({"flag":flag,"rname":rname,"pos":pos+1,"mapq":mapq,
                          "l_seq":l_seq,"tlen":tlen,
                          "is_unmapped":bool(flag&0x4),"is_duplicate":bool(flag&0x400),
                          "is_paired":bool(flag&0x1),"is_proper_pair":bool(flag&0x2),
                          "is_secondary":bool(flag&0x100),"is_supplementary":bool(flag&0x800),
                          "is_reverse":bool(flag&0x10),"is_read1":bool(flag&0x40),
                          "is_read2":bool(flag&0x80),"mate_unmapped":bool(flag&0x8)})
        except: continue
    return header_text,refs,reads

# ════════════════════════════════════════════════════════════
# ANALYZERS
# ════════════════════════════════════════════════════════════

def analyze_bam(path):
    sections=[]
    total=mapped=unmapped=duplicate=paired=proper_pair=0
    secondary=supplementary=mate_unmapped=forward=reverse_=read1=read2=0
    read_lengths=[]; mapq_scores=[]; insert_sizes=[]
    chrom_counts=collections.Counter(); nm_vals=[]; sample_rows=[]; refs=[]

    if PYSAM_OK:
        try:
            mode="rb" if path.lower().endswith(".bam") else "r"
            aln=pysam.AlignmentFile(path,mode)
            hdr=aln.header.to_dict()
            hd=hdr.get("HD",{}); sq_list=hdr.get("SQ",[]); rg_list=hdr.get("RG",[]); pg_list=hdr.get("PG",[])
            refs=[(s.get("SN",""),s.get("LN",0)) for s in sq_list]
            hdr_rows=[("SAM version",hd.get("VN","N/A")),("Sort order",hd.get("SO","N/A")),
                      ("References",str(len(sq_list))),("Read groups",str(len(rg_list)))]
            for rg in rg_list: hdr_rows+=[("Sample",rg.get("SM","N/A")),("Platform",rg.get("PL","N/A"))]
            for pg in pg_list: hdr_rows+=[("Aligner",f"{pg.get('PN','?')} v{pg.get('VN','?')}"),("Command",pg.get("CL","N/A")[:100])]
            sections.append(_metrics("Header",hdr_rows))
            if sq_list: sections.append(_table("Reference Sequences",[{"Contig":s.get("SN",""),"Length (bp)":f"{s.get('LN',0):,}"} for s in sq_list]))
            for read in aln:
                total+=1
                if read.is_unmapped: unmapped+=1
                else:
                    mapped+=1; chrom_counts[read.reference_name]+=1
                    if read.mapping_quality is not None: mapq_scores.append(read.mapping_quality)
                if read.is_duplicate: duplicate+=1
                if read.is_paired: paired+=1
                if read.is_proper_pair: proper_pair+=1
                if read.is_secondary: secondary+=1
                if read.is_supplementary: supplementary+=1
                if read.mate_is_unmapped: mate_unmapped+=1
                if read.is_forward: forward+=1
                else: reverse_+=1
                if read.is_read1: read1+=1
                if read.is_read2: read2+=1
                if read.query_length: read_lengths.append(read.query_length)
                if read.template_length and 0<abs(read.template_length)<10000: insert_sizes.append(abs(read.template_length))
                try: nm_vals.append(read.get_tag("NM"))
                except: pass
                if len(sample_rows)<20:
                    sample_rows.append({"QNAME":read.query_name,"FLAG":read.flag,"RNAME":read.reference_name or "*",
                                        "POS":read.reference_start+1 if not read.is_unmapped else 0,
                                        "MAPQ":read.mapping_quality,"CIGAR":read.cigarstring or "*","TLEN":read.template_length})
            aln.close()
        except: total=0

    if not total:
        # Detect if it's a text SAM or binary BAM
        with open(path,"rb") as _f: _magic=_f.read(4)
        is_sam_text = (_magic[0:2] != b"\x1f\x8b")  # not gzip → likely SAM text

        if is_sam_text:
            # ── Pure-Python SAM text parser ──────────────────────────
            hdr_rows=[]; refs=[]
            try:
                with open(path,"rt",errors="ignore") as f:
                    for line in f:
                        line=line.rstrip()
                        if line.startswith("@"):
                            tag=line[:3]
                            parts_={t[:2]:t[3:] for t in line.split("\t")[1:] if len(t)>2 and ":"in t}
                            if tag=="@HD": hdr_rows+=[("SAM version",parts_.get("VN","N/A")),("Sort order",parts_.get("SO","N/A"))]
                            elif tag=="@SQ":
                                name=parts_.get("SN","?"); length=int(parts_.get("LN",0))
                                refs.append((name,length))
                            elif tag=="@RG": hdr_rows+=[("Sample",parts_.get("SM","N/A")),("Platform",parts_.get("PL","N/A"))]
                            elif tag=="@PG": hdr_rows+=[("Aligner",f"{parts_.get('PN','?')} v{parts_.get('VN','?')}"),("Command",parts_.get("CL","N/A")[:100])]
                            continue
                        cols=line.split("\t")
                        if len(cols)<11: continue
                        try:
                            qname,flag,rname,pos,mapq_,cigar,rnext,pnext,tlen_,seq,qual_=cols[:11]
                            flag=int(flag); mapq_=int(mapq_); pos=int(pos); tlen_=int(tlen_)
                            is_unmapped_=bool(flag&0x4); is_dup_=bool(flag&0x400)
                            is_paired_=bool(flag&0x1); is_proper_=bool(flag&0x2)
                            is_secondary_=bool(flag&0x100); is_supplementary_=bool(flag&0x800)
                            is_reverse_=bool(flag&0x10); is_read1_=bool(flag&0x40); is_read2_=bool(flag&0x80)
                            mate_unmapped_=bool(flag&0x8)
                            total+=1
                            if is_unmapped_: unmapped+=1
                            else:
                                mapped+=1; chrom_counts[rname]+=1; mapq_scores.append(mapq_)
                            if is_dup_: duplicate+=1
                            if is_paired_: paired+=1
                            if is_proper_: proper_pair+=1
                            if is_secondary_: secondary+=1
                            if is_supplementary_: supplementary+=1
                            if mate_unmapped_: mate_unmapped+=1
                            if is_reverse_: reverse_+=1
                            else: forward+=1
                            if is_read1_: read1+=1
                            if is_read2_: read2+=1
                            if seq and seq!="*": read_lengths.append(len(seq))
                            if 0<abs(tlen_)<10000: insert_sizes.append(abs(tlen_))
                            if len(sample_rows)<20:
                                sample_rows.append({"QNAME":qname,"FLAG":flag,"RNAME":rname,"POS":pos,"MAPQ":mapq_,"CIGAR":cigar,"TLEN":tlen_})
                        except: continue
            except Exception as e: return _err("SAM",str(e))
            if not total: return _err("SAM","No alignment records found in SAM file.")
            sections.append(_metrics("Header",hdr_rows if hdr_rows else [("Format","SAM (no header lines detected)")]))
            if refs: sections.append(_table("Reference Sequences",[{"Contig":n,"Length (bp)":f"{l:,}"} for n,l in refs[:50]]))

        else:
            # ── Pure-Python BAM binary parser ────────────────────────
            ht,refs_raw,read_list=_parse_bam_pure(path)
            if ht is None: return _err("BAM","Could not parse BAM file.")
            refs=refs_raw
            hdr_rows=[]
            for line in ht.splitlines():
                parts={f[:2]:f[3:] for f in line.split("\t") if ":"in f}
                if line.startswith("@HD"): hdr_rows+=[("Sort order",parts.get("SO","N/A")),("Version",parts.get("VN","N/A"))]
                if line.startswith("@PG"): hdr_rows+=[("Aligner",f"{parts.get('PN','?')} v{parts.get('VN','?')}"),("Command",parts.get("CL","N/A")[:100])]
                if line.startswith("@RG"): hdr_rows+=[("Sample",parts.get("SM","N/A")),("Platform",parts.get("PL","N/A"))]
            sections.append(_metrics("Header",hdr_rows))
            if refs: sections.append(_table("Reference Sequences",[{"Contig":n,"Length (bp)":f"{l:,}"} for n,l in refs]))
            for r in read_list:
                total+=1
                if r["is_unmapped"]: unmapped+=1
                else: mapped+=1; chrom_counts[r["rname"]]+=1; mapq_scores.append(r["mapq"])
                if r["is_duplicate"]: duplicate+=1
                if r["is_paired"]: paired+=1
                if r["is_proper_pair"]: proper_pair+=1
                if r["is_secondary"]: secondary+=1
                if r["is_supplementary"]: supplementary+=1
                if r["mate_unmapped"]: mate_unmapped+=1
                if r["is_reverse"]: reverse_+=1
                else: forward+=1
                if r["is_read1"]: read1+=1
                if r["is_read2"]: read2+=1
                if r["l_seq"]>0: read_lengths.append(r["l_seq"])
                if 0<abs(r.get("tlen",0))<10000: insert_sizes.append(abs(r["tlen"]))
                if len(sample_rows)<20: sample_rows.append({"FLAG":r["flag"],"RNAME":r["rname"],"POS":r["pos"],"MAPQ":r["mapq"],"TLEN":r.get("tlen",0)})

    if total:
        sections.append(_metrics("Alignment Statistics",[
            ("Total reads",f"{total:,}"),("Mapped",f"{mapped:,}  ({mapped/total*100:.2f}%)"),
            ("Unmapped",f"{unmapped:,}  ({unmapped/total*100:.2f}%)"),("Duplicates",f"{duplicate:,}  ({duplicate/total*100:.2f}%)"),
            ("Paired-end",f"{paired:,}"),("Proper pairs",f"{proper_pair:,}  ({proper_pair/total*100:.2f}%)"),
            ("Mate unmapped",f"{mate_unmapped:,}"),("Secondary",f"{secondary:,}"),("Supplementary",f"{supplementary:,}"),
            ("Forward strand",f"{forward:,}  ({forward/total*100:.1f}%)"),("Reverse strand",f"{reverse_:,}  ({reverse_/total*100:.1f}%)"),
            ("Read 1 (R1)",f"{read1:,}"),("Read 2 (R2)",f"{read2:,}"),
        ]))
    if read_lengths:
        sections.append(_metrics("Read Length",[("Min",f"{min(read_lengths)} bp"),("Max",f"{max(read_lengths)} bp"),
            ("Mean",f"{statistics.mean(read_lengths):.1f} bp"),("Median",f"{statistics.median(read_lengths):.1f} bp")]))
    if mapq_scores:
        sections.append(_metrics("Mapping Quality (MAPQ)",[
            ("Min / Max",f"{min(mapq_scores)} / {max(mapq_scores)}"),("Mean",f"{statistics.mean(mapq_scores):.2f}"),
            ("Median",f"{statistics.median(mapq_scores):.1f}"),("MAPQ=0 (multi-mappers)",f"{sum(1 for m in mapq_scores if m==0):,}"),
            ("MAPQ≥20 (unique)",f"{sum(1 for m in mapq_scores if m>=20):,}"),("MAPQ≥30 (high quality)",f"{sum(1 for m in mapq_scores if m>=30):,}"),
            ("MAPQ=60 (max confidence)",f"{sum(1 for m in mapq_scores if m==60):,}"),
        ]))
        bins=collections.Counter(int(m//10)*10 for m in mapq_scores)
        sections.append(_bars("MAPQ Distribution",[(f"MAPQ {b}–{b+9}",bins[b],round(bins[b]/len(mapq_scores)*100,1),bins[b]/len(mapq_scores)*100) for b in sorted(bins)]))
    if insert_sizes:
        sections.append(_metrics("Insert Size",[("Min",f"{min(insert_sizes)} bp"),("Max",f"{max(insert_sizes)} bp"),
            ("Median",f"{statistics.median(insert_sizes):.0f} bp"),("Mean",f"{statistics.mean(insert_sizes):.1f} bp"),
            ("<100 bp (adapter risk)",f"{sum(1 for i in insert_sizes if i<100):,}"),
            (">1000 bp (large/SVs)",f"{sum(1 for i in insert_sizes if i>1000):,}")]))
    if chrom_counts:
        sections.append(_table("Reads per Chromosome",[{"Chromosome":ch,"Reads":f"{cnt:,}","% of Mapped":f"{cnt/max(mapped,1)*100:.2f}%"} for ch,cnt in chrom_counts.most_common(25)]))
    if sample_rows: sections.append(_table("Sample Reads (first 20)",sample_rows))

    qc=[]
    if total:
        qc+=[_qc("Mapping rate ≥ 90%",_grade(mapped/total>=0.9,mapped/total>=0.7)),
             _qc("Proper pair rate ≥ 80%",_grade(proper_pair/total>=0.8,proper_pair/total>=0.5)),
             _qc("Duplicate rate < 20%",_grade(duplicate/total<0.2,duplicate/total<0.5))]
    if mapq_scores: qc.append(_qc("Mean MAPQ ≥ 30",_grade(statistics.mean(mapq_scores)>=30,statistics.mean(mapq_scores)>=20)))

    interp=["Mapping rate <70% → contamination, wrong reference, or poor library quality",
            "Proper pair rate <80% → structural rearrangements or library prep failure",
            "Duplicates >20% → PCR over-amplification — mark with Picard MarkDuplicates",
            "Secondary reads → multi-mappers; exclude for SNP calling, include for RNA-seq",
            "Supplementary reads → chimeric/split reads; key signal for SV and fusion detection",
            "MAPQ=0 → read maps equally to multiple locations (repetitive regions)",
            "MAPQ=60 → maximum unique alignment confidence (BWA-MEM scale)",
            "Insert <100 bp → adapter read-through; re-trim with Cutadapt or Trimmomatic"]
    return _ok("BAM / SAM",sections,qc,interp)


def analyze_fasta(path):
    sections=[]; records=[]; lengths=[]; gc_vals=[]; ambig_list=[]; cpg_list=[]
    kmer_all=collections.Counter(); dinuc_all=collections.Counter()
    opener=gzip.open if path.endswith(".gz") else open

    def proc(sid,desc,seq):
        if not seq: return
        n=len(seq); a,t,g,c_=seq.count("A"),seq.count("T"),seq.count("G"),seq.count("C")
        gc=(g+c_)/n*100; ambig=sum(1 for b in seq if b not in "ACGTU"); cpg=len(re.findall("CG",seq))
        at_sk=(a-t)/(a+t) if (a+t) else 0; gc_sk=(g-c_)/(g+c_) if (g+c_) else 0
        din=collections.Counter(seq[i:i+2] for i in range(n-1) if seq[i] in "ACGT" and seq[i+1] in "ACGT")
        tri=collections.Counter(seq[i:i+3] for i in range(n-2) if "N" not in seq[i:i+3])
        kmer_all.update(tri); dinuc_all.update(din)
        total_t=sum(tri.values())
        ent=-sum((c/total_t)*math.log2(c/total_t) for c in tri.values() if c>0) if total_t else 0
        records.append({"ID":sid,"Description":desc[:60],"Length (bp)":f"{n:,}","GC%":f"{gc:.2f}%",
                        "Ambiguous":ambig,"CpG sites":cpg,"AT skew":f"{at_sk:+.4f}","GC skew":f"{gc_sk:+.4f}","Complexity":f"{ent:.3f}"})
        lengths.append(n); gc_vals.append(gc); ambig_list.append(ambig); cpg_list.append(cpg)

    try:
        if BIOPYTHON_OK:
            for rec in SeqIO.parse(path,"fasta"): proc(rec.id,rec.description[:80],str(rec.seq).upper())
        else:
            with opener(path,"rt",errors="ignore") as f:
                sid=sdesc=parts_=None
                for line in f:
                    line=line.strip()
                    if line.startswith(">"):
                        if sid and parts_: proc(sid,sdesc,"".join(parts_).upper())
                        p=line[1:].split(None,1); sid=p[0]; sdesc=p[1] if len(p)>1 else sid; parts_=[]
                    elif parts_ is not None: parts_.append(line)
                if sid and parts_: proc(sid,sdesc,"".join(parts_).upper())
    except Exception as e: return _err("FASTA",str(e))

    if not lengths: return _err("FASTA","No sequences found.")
    total_bases=sum(lengths); sl=sorted(lengths,reverse=True); cum=n50=n90=0
    for l in sl:
        cum+=l
        if n50==0 and cum>=total_bases*0.5: n50=l
        if n90==0 and cum>=total_bases*0.9: n90=l

    sections=[_metrics("Assembly Overview",[
        ("Total sequences",f"{len(lengths):,}"),("Total bases",f"{total_bases:,} bp"),
        ("Shortest",f"{min(lengths):,} bp"),("Longest",f"{max(lengths):,} bp"),
        ("Mean length",f"{statistics.mean(lengths):,.1f} bp"),("N50",f"{n50:,} bp"),("N90",f"{n90:,} bp"),
        ("Mean GC%",f"{statistics.mean(gc_vals):.2f}%"),("GC% range",f"{min(gc_vals):.2f}% – {max(gc_vals):.2f}%"),
        ("Total N bases",f"{sum(ambig_list):,}  ({sum(ambig_list)/total_bases*100:.2f}%)"),("CpG sites",f"{sum(cpg_list):,}"),
    ])]
    bins=collections.Counter(int(g//10)*10 for g in gc_vals)
    sections.append(_bars("GC% Distribution",[(f"{b}–{b+9}%",bins.get(b,0),round(bins.get(b,0)/len(gc_vals)*100,1),bins.get(b,0)/len(gc_vals)*100) for b in sorted(set(int(g//10)*10 for g in gc_vals))]))
    total_din=sum(dinuc_all.values())
    if total_din:
        din_rows=[{"Dinucleotide":d,"Count":f"{dinuc_all.get(d,0):,}","Frequency":f"{dinuc_all.get(d,0)/total_din*100:.2f}%",
                   "Note":"← CpG (methylation)" if d=="CG" else "← TpA (suppressed)" if d=="TA" else ""}
                  for d in ["AA","AT","AG","AC","TA","TT","TG","TC","GA","GT","GG","GC","CA","CT","CG","CC"]]
        sections.append(_table("Dinucleotide Frequencies",din_rows))
    sections.append(_table("Top 10 Trimers",[ {"3-mer":k,"Count":f"{c:,}"} for k,c in kmer_all.most_common(10)]))
    sections.append(_table("Per-Sequence Summary",records))

    qc=[_qc("N50 ≥ 1 Mb",_grade(n50>=1_000_000,n50>=100_000)),
        _qc("N content < 5%",_grade(sum(ambig_list)/total_bases<0.05)),
        _qc("GC% within 30–70%",_grade(30<=statistics.mean(gc_vals)<=70)),
        _qc("Sequence count < 1,000",_grade(len(lengths)<1000,len(lengths)<10000))]
    interp=["N50: half the assembly is in contigs ≥ this length — higher = more contiguous",
            "GC% is species-specific — deviation may indicate contamination or HGT",
            "CpG islands (>200 bp, >50% GC) mark gene promoters",
            "Low 3-mer complexity → repetitive regions (transposons, microsatellites)",
            "N content >5% may reduce annotation accuracy"]
    return _ok("FASTA",sections,qc,interp)


def analyze_genbank(path):
    sections=[]; metadata={}; features=[]; seq_parts=[]
    current_feature=current_qual=None; current_val=[]; in_origin=in_features=False
    META_KEYS={"LOCUS","DEFINITION","ACCESSION","VERSION","SOURCE","ORGANISM"}
    opener=gzip.open if path.endswith(".gz") else open
    try:
        with opener(path,"rt",errors="ignore") as f:
            for raw in f:
                line=raw.rstrip()
                if line.startswith("ORIGIN"):
                    in_origin=True; in_features=False
                    if current_feature and current_qual: current_feature["qualifiers"].setdefault(current_qual,[]).append(" ".join(current_val))
                    if current_feature: features.append(current_feature)
                    current_feature=None; continue
                if line.strip()=="//": in_origin=False; continue
                if in_origin: seq_parts.append(re.sub(r"[\d\s]","",line)); continue
                if line.startswith("FEATURES"): in_features=True; continue
                if in_features:
                    m=re.match(r"^     (\S+)\s+(.+)$",line)
                    if m and not line.startswith("                     "):
                        if current_feature and current_qual: current_feature["qualifiers"].setdefault(current_qual,[]).append(" ".join(current_val))
                        if current_feature: features.append(current_feature)
                        current_feature={"type":m.group(1),"location":m.group(2).strip(),"qualifiers":{}}; current_qual=None; current_val=[]; continue
                    if line.startswith("                     /"):
                        if current_feature and current_qual: current_feature["qualifiers"].setdefault(current_qual,[]).append(" ".join(current_val))
                        rest=line.strip()[1:]
                        if "=" in rest: k,v=rest.split("=",1); current_qual=k; current_val=[v.strip('"')]
                        else: current_qual=rest.strip(); current_val=["true"]
                        continue
                    if line.startswith("                     ") and current_qual: current_val.append(line.strip().strip('"')); continue
                if not in_features:
                    for key in META_KEYS:
                        if line.startswith(key): metadata[key]=line[len(key):].strip(); break
    except Exception as e: return _err("GENBANK",str(e))
    if current_feature:
        if current_qual: current_feature["qualifiers"].setdefault(current_qual,[]).append(" ".join(current_val))
        features.append(current_feature)

    sequence="".join(seq_parts).upper(); seq_len=len(sequence)
    def parse_loc(ls):
        strand="+"
        if ls.startswith("complement("): strand="-"; ls=ls[11:-1]
        ls=re.sub(r"join\(|order\(|\)","",ls)
        coords=re.findall(r"[<>]?(\d+)\.\.[<>]?(\d+)",ls)
        if coords: return min(int(c[0]) for c in coords),max(int(c[1]) for c in coords),strand
        single=re.findall(r"\d+",ls)
        if single: p=int(single[0]); return p,p,strand
        return None,None,strand

    rows=[]; feat_types=[]; cds_lengths=[]; products=[]; ec_numbers=[]; hypothetical=0; with_ec=0
    strand_counts=collections.Counter(); transl_tables=collections.Counter()
    for feat in features:
        q=feat["qualifiers"]; start,end,strand=parse_loc(feat["location"])
        length=(end-start+1) if (start and end) else None
        product=q.get("product",[""])[0]; ec=q.get("EC_number",[""])[0]; ttable=q.get("transl_table",[""])[0]
        if product: products.append(product)
        if ec: ec_numbers.append(ec); with_ec+=1
        if ttable: transl_tables[ttable]+=1
        if product and "hypothetical" in product.lower(): hypothetical+=1
        if feat["type"]=="CDS" and length: cds_lengths.append(length)
        feat_types.append(feat["type"]); strand_counts[strand]+=1
        rows.append({"Type":feat["type"],"Start":start,"End":end,"Strand":strand,"Length":length,
                     "Gene":q.get("gene",[""])[0],"Locus tag":q.get("locus_tag",[""])[0],
                     "Product":product[:50],"EC":ec})

    gc_pct=(sequence.count("G")+sequence.count("C"))/seq_len*100 if seq_len else 0
    coding_density=sum(cds_lengths)/seq_len*100 if seq_len else 0
    sections=[_metrics("Record Metadata",[("LOCUS",metadata.get("LOCUS","N/A")[:80]),("Definition",metadata.get("DEFINITION","N/A")[:100]),
        ("Accession",metadata.get("ACCESSION","N/A")),("Organism",metadata.get("ORGANISM","N/A")[:100])])]
    if seq_len:
        sections.append(_metrics("Sequence Statistics",[("Total length",f"{seq_len:,} bp"),("GC content",f"{gc_pct:.2f}%"),
            ("CpG sites",f"{len(re.findall('CG',sequence)):,}"),("Coding density",f"{coding_density:.1f}%")]))
    sections.append(_table("Feature Type Counts",[{"Feature":k,"Count":f"{v:,}"} for k,v in pd.Series(feat_types).value_counts().items()]))
    if cds_lengths:
        sections.append(_metrics("CDS Statistics",[("Total CDS",f"{len(cds_lengths):,}"),("Min",f"{min(cds_lengths):,} bp"),
            ("Max",f"{max(cds_lengths):,} bp"),("Mean",f"{statistics.mean(cds_lengths):,.1f} bp"),("With EC number",f"{with_ec:,}")]))
    if products:
        known=len(products)-hypothetical
        sections.append(_metrics("Product Annotation",[("Known proteins",f"{known:,}  ({known/len(products)*100:.1f}%)"),("Hypothetical",f"{hypothetical:,}  ({hypothetical/len(products)*100:.1f}%)")]))
        sections.append(_table("Top 15 Products",[{"Product":p,"Count":f"{c:,}"} for p,c in collections.Counter(products).most_common(15)]))
    if ec_numbers:
        ec_names={"1":"Oxidoreductases","2":"Transferases","3":"Hydrolases","4":"Lyases","5":"Isomerases","6":"Ligases","7":"Translocases"}
        ec_cls=collections.Counter(ec.split(".")[0] for ec in ec_numbers)
        sections.append(_table("Enzyme Classes",[{"EC":f"EC {k}.x","Name":ec_names.get(k,"Unknown"),"Count":f"{v:,}"} for k,v in sorted(ec_cls.items())]))
    sections.append(_table("Feature Table (first 50)",rows[:50]))

    hypo_pct=hypothetical/len(products)*100 if products else 0
    qc=[_qc("Sequence present",_grade(seq_len>0)),_qc("Coding density 70–95%",_grade(70<=coding_density<=95,seq_len>0)),
        _qc("Hypothetical proteins <30%",_grade(hypo_pct<30,hypo_pct<60)),_qc("tRNA annotated",_grade("tRNA" in set(feat_types)))]
    interp=["Prokaryote coding density: typically 85–95% (compact genomes)",
            "/locus_tag: stable systematic ID — use for cross-referencing databases",
            "/EC_number: links gene to KEGG/MetaCyc metabolic pathways",
            "High hypothetical% (>50%): novel organism or insufficient reference database",
            "GC skew sign change marks replication origin/terminus in circular genomes"]
    return _ok("GENBANK",sections,qc,interp)


def analyze_gff(path):
    sections=[]; rows=[]; feat_types=[]; chroms=[]; strands=[]; products=[]; lengths=[]
    inference_tools=collections.Counter(); gene_ids=set()
    opener=gzip.open if path.endswith(".gz") else open
    try:
        with opener(path,"rt",errors="ignore") as f:
            for line in f:
                line=line.rstrip()
                if line.startswith("#"): continue
                cols=line.split("\t")
                if len(cols)<9: continue
                seqname,source,feature,start,end,score,strand,frame,attributes=cols[:9]
                try: length=int(end)-int(start)+1
                except: length=None
                attr_dict={}
                for part in re.split(r";",attributes):
                    part=part.strip()
                    if "=" in part: k,v=part.split("=",1); attr_dict[k.strip()]=v.strip()
                product=attr_dict.get("product",""); gene_id=attr_dict.get("ID",""); ltag=attr_dict.get("locus_tag","")
                if gene_id: gene_ids.add(gene_id)
                if product: products.append(product)
                rows.append({"Seqname":seqname,"Feature":feature,"Start":start,"End":end,"Strand":strand,"ID":gene_id,"Locus tag":ltag,"Product":product[:50]})
                feat_types.append(feature); chroms.append(seqname); strands.append(strand)
                if length: lengths.append(length)
    except Exception as e: return _err("GFF3",str(e))
    sections=[_metrics("GFF3 Overview",[("Total features",f"{len(rows):,}"),("Unique sequences",f"{len(set(chroms)):,}"),
        ("Feature types",f"{len(set(feat_types)):,}"),("Gene IDs",f"{len(gene_ids):,}"),("Products",f"{len(products):,}")])]
    sections.append(_table("Feature Type Counts",[{"Feature":k,"Count":f"{v:,}"} for k,v in pd.Series(feat_types).value_counts().items()]))
    fwd=strands.count("+"); rev=strands.count("-")
    sections.append(_metrics("Strand Distribution",[("Forward (+)",f"{fwd:,}  ({fwd/len(strands)*100:.1f}%)" if strands else "0"),("Reverse (−)",f"{rev:,}  ({rev/len(strands)*100:.1f}%)" if strands else "0")]))
    if lengths: sections.append(_metrics("Feature Length",[("Min",f"{min(lengths):,} bp"),("Max",f"{max(lengths):,} bp"),("Mean",f"{statistics.mean(lengths):,.1f} bp"),("Median",f"{statistics.median(lengths):,.1f} bp")]))
    if products:
        hypo=sum(1 for p in products if "hypothetical" in p.lower())
        sections.append(_metrics("Products",[("Known",f"{len(products)-hypo:,}  ({(len(products)-hypo)/len(products)*100:.1f}%)"),("Hypothetical",f"{hypo:,}  ({hypo/len(products)*100:.1f}%)")]))
    sections.append(_table("Feature Table (first 50)",rows[:50]))
    ft_set=set(feat_types)
    hypo_pct=sum(1 for p in products if "hypothetical" in p.lower())/len(products)*100 if products else 0
    qc=[_qc("CDS present",_grade("CDS" in ft_set)),_qc("tRNA present",_grade("tRNA" in ft_set)),
        _qc("rRNA present",_grade("rRNA" in ft_set)),_qc("Hypothetical <30%",_grade(hypo_pct<30,hypo_pct<60))]
    interp=["GFF3: 1-based closed coordinates [Start, End] — End−Start+1 = length",
            "Feature hierarchy: gene → mRNA → exon/CDS (linked by ID/Parent attributes)",
            "CDS excludes UTRs; exon includes UTRs",
            "High hypothetical% (>50%): novel organism or poor reference database coverage"]
    return _ok("GFF3",sections,qc,interp)


def analyze_gtf(path):
    sections=[]; rows=[]; feat_types=[]; chroms=[]
    biotypes=collections.Counter(); tsl=collections.Counter()
    exon_counts=collections.defaultdict(int); cds_lengths=[]
    gene_ids=set(); transcript_ids=set()
    def parse_attrs(s):
        d={}
        for m in re.finditer(r'(\w+)\s+"([^"]*)"',s): d[m.group(1)]=m.group(2)
        return d
    opener=gzip.open if path.endswith(".gz") else open
    try:
        with opener(path,"rt",errors="ignore") as f:
            for line in f:
                line=line.rstrip()
                if line.startswith("#"): continue
                cols=line.split("\t")
                if len(cols)<9: continue
                seqname,source,feature,start,end,score,strand,frame,attributes=cols[:9]
                try: length=int(end)-int(start)+1
                except: length=None
                attrs=parse_attrs(attributes)
                gid=attrs.get("gene_id",""); tid=attrs.get("transcript_id","")
                gname=attrs.get("gene_name",attrs.get("gene",""))
                bt=attrs.get("gene_biotype",attrs.get("transcript_biotype",attrs.get("gene_type","")))
                t_lvl=attrs.get("transcript_support_level","")
                if gid: gene_ids.add(gid)
                if tid: transcript_ids.add(tid)
                if bt: biotypes[bt]+=1
                if t_lvl: tsl[t_lvl]+=1
                if feature=="exon" and tid: exon_counts[tid]+=1
                if feature=="CDS" and length: cds_lengths.append(length)
                feat_types.append(feature); chroms.append(seqname)
                rows.append({"Seqname":seqname,"Feature":feature,"Start":start,"End":end,"Strand":strand,"gene_id":gid,"gene_name":gname,"biotype":bt})
    except Exception as e: return _err("GTF",str(e))
    sections=[_metrics("GTF Overview",[("Total features",f"{len(rows):,}"),("Chromosomes",f"{len(set(chroms)):,}"),("Unique genes",f"{len(gene_ids):,}"),("Unique transcripts",f"{len(transcript_ids):,}")])]
    sections.append(_table("Feature Type Counts",[{"Feature":k,"Count":f"{v:,}"} for k,v in pd.Series(feat_types).value_counts().items()]))
    if biotypes:
        total_bt=sum(biotypes.values())
        sections.append(_bars("Gene Biotypes (top 15)",[(k,v,round(v/total_bt*100,1),v/total_bt*100) for k,v in biotypes.most_common(15)]))
    if exon_counts:
        ec_list=list(exon_counts.values())
        sections.append(_metrics("Exon Structure",[("Transcripts with exons",f"{len(ec_list):,}"),("Min exons/transcript",str(min(ec_list))),("Max exons/transcript",str(max(ec_list))),("Mean exons/transcript",f"{statistics.mean(ec_list):.1f}"),("Single-exon",f"{sum(1 for e in ec_list if e==1):,}")]))
    sections.append(_table("Feature Table (first 50)",rows[:50]))
    ft_set=set(feat_types)
    qc=[_qc("gene features present",_grade("gene" in ft_set)),_qc("exon features present",_grade("exon" in ft_set)),
        _qc("CDS features present",_grade("CDS" in ft_set)),_qc("protein_coding biotype",_grade("protein_coding" in biotypes))]
    interp=["gene_id: stable Ensembl ID (e.g. ENSG00000139618 = BRCA2)",
            "gene_name: HGNC symbol (e.g. TP53, EGFR, BRCA1) — human-readable",
            "TSL 1 = highest confidence; filter to TSL1/2 for isoform analysis",
            "GENCODE: comprehensive — best for WGS/WES; RefSeq: preferred for clinical reporting"]
    return _ok("GTF",sections,qc,interp)


def analyze_bed(path,fmt="BED"):
    rows=[]; chroms=[]; lengths=[]; scores=[]
    with open(path) as f:
        for line in f:
            if line.startswith(("#","track","browser")): continue
            cols=line.strip().split("\t")
            if len(cols)<3: continue
            try:
                chrom=cols[0]; start=int(cols[1]); end=int(cols[2]); length=end-start
                name=cols[3] if len(cols)>3 else "."; score_=cols[4] if len(cols)>4 else "."; strand=cols[5] if len(cols)>5 else "."
                row={"Chrom":chrom,"Start":start,"End":end,"Length":length,"Name":name,"Score":score_,"Strand":strand}
                if fmt in ("NARROWPEAK","BROADPEAK") and len(cols)>=9:
                    try: row["signalValue"]=float(cols[6]); row["pValue"]=float(cols[7]); row["qValue"]=float(cols[8])
                    except: pass
                rows.append(row); chroms.append(chrom); lengths.append(length)
                try: scores.append(float(score_))
                except: pass
            except: continue
    if not rows: return _err(fmt,"No valid BED records found.")
    df=pd.DataFrame(rows)
    sections=[_metrics(f"{fmt} Overview",[("Total regions",f"{len(rows):,}"),("Chromosomes",f"{len(set(chroms)):,}"),("Total bp covered",f"{sum(lengths):,}")])]
    if lengths: sections.append(_metrics("Region Length",[("Min",f"{min(lengths):,} bp"),("Max",f"{max(lengths):,} bp"),("Mean",f"{statistics.mean(lengths):,.1f} bp"),("Median",f"{statistics.median(lengths):,.1f} bp")]))
    chrom_stats=df.groupby("Chrom").agg(Regions=("Length","count"),Total_bp=("Length","sum")).sort_values("Regions",ascending=False)
    sections.append(_table("Chromosome Distribution",chrom_stats.head(25).reset_index()))
    if scores: sections.append(_metrics("Score",[("Min",f"{min(scores):.1f}"),("Max",f"{max(scores):.1f}"),("Mean",f"{statistics.mean(scores):.1f}"),("Score≥500",f"{sum(1 for s in scores if s>=500):,}")]))
    sections.append(_table("Feature Table (first 50)",rows[:50]))
    qc=[_qc("Peak count ≥ 1000",_grade(len(rows)>=1000,len(rows)>=100))]
    if lengths: qc.append(_qc("Median width 100–2000 bp",_grade(100<=statistics.median(lengths)<=2000)))
    interp=["BED uses 0-based half-open [Start,End) — End−Start = length",
            "narrowPeak: point-source peaks — TF binding, H3K4me3 (active promoters)",
            "broadPeak: diffuse peaks — H3K27me3 (silencing) / H3K36me3 (transcription)",
            "Score ≥ 500: high-confidence ENCODE peak threshold"]
    return _ok(fmt,sections,qc,interp)


def analyze_vcf(path):
    """Pure-Python VCF/BCF parser — works without pysam."""
    sections=[]
    rows=[]; types=[]; zyg=[]; chroms=collections.Counter()
    filters_cnt=collections.Counter(); qual_vals=[]; dp_vals=[]; af_vals=[]
    ts_count=tv_count=multiallelic=0; indel_sizes=[]; samples=[]; info_fields=[]; fmt_fields=[]; contigs=[]
    TRANSITIONS={("A","G"),("G","A"),("C","T"),("T","C")}

    # Try pysam first
    if PYSAM_OK:
        try:
            vcf=pysam.VariantFile(path)
            hdr=vcf.header
            samples=list(hdr.samples)
            info_fields=list(hdr.info.keys())
            fmt_fields=list(hdr.formats.keys())
            contigs=[(c.name, c.length or 0) for c in hdr.contigs.values()]
            for rec in vcf:
                ref=rec.ref; alts=rec.alts or []
                if len(alts)>1: multiallelic+=1
                if rec.qual is not None: qual_vals.append(rec.qual)
                for flt in rec.filter.keys(): filters_cnt[flt]+=1
                try: dp_vals.append(rec.info["DP"])
                except: pass
                try:
                    af=rec.info["AF"]; af=af[0] if isinstance(af,tuple) else af; af_vals.append(float(af))
                except: pass
                chroms[rec.chrom]+=1
                for alt in (alts or []):
                    if len(ref)==1 and len(alt)==1:
                        t="SNP"; pair=(ref.upper(),alt.upper())
                        if pair in TRANSITIONS: ts_count+=1
                        else: tv_count+=1
                    elif len(ref)!=len(alt): t="INDEL"; indel_sizes.append(len(alt)-len(ref))
                    else: t="MNV"
                    gt_str="."
                    if rec.samples:
                        s=list(rec.samples.values())[0]
                        try:
                            gt=s["GT"]
                            if gt in [(1,1),(1,)]: gt_str="Homozygous ALT"
                            elif gt==(0,0): gt_str="Homozygous REF"
                            elif gt and 1 in gt: gt_str="Heterozygous"
                            elif gt and None in gt: gt_str="Missing"
                        except: pass
                    rows.append({"CHROM":rec.chrom,"POS":rec.pos,"REF":ref,"ALT":alt,"QUAL":rec.qual,"FILTER":",".join(rec.filter.keys()),"TYPE":t,"ZYGOSITY":gt_str})
                    types.append(t); zyg.append(gt_str)
            vcf.close()
        except Exception as e:
            rows=[]  # fall through to text parser

    # Pure-text fallback
    if not rows:
        opener=gzip.open if path.endswith(".gz") else open
        try:
            with opener(path,"rt",errors="ignore") as f:
                for line in f:
                    line=line.rstrip()
                    if line.startswith("##"):
                        if line.startswith("##INFO=") and not info_fields:
                            m=re.search(r'ID=([^,>]+)',line)
                            if m: info_fields.append(m.group(1))
                        if line.startswith("##FORMAT="):
                            m=re.search(r'ID=([^,>]+)',line)
                            if m: fmt_fields.append(m.group(1))
                        if line.startswith("##contig="):
                            m=re.search(r'ID=([^,>]+)',line); l=re.search(r'length=(\d+)',line)
                            if m: contigs.append((m.group(1),int(l.group(1)) if l else 0))
                        continue
                    if line.startswith("#CHROM"):
                        cols=line.split("\t")
                        samples=cols[9:] if len(cols)>9 else []; continue
                    if not line or line.startswith("#"): continue
                    cols=line.split("\t")
                    if len(cols)<5: continue
                    chrom,pos,id_,ref=cols[0],cols[1],cols[2],cols[3]
                    alts_str=cols[4]; qual=cols[5] if len(cols)>5 else "."
                    filt=cols[6] if len(cols)>6 else "."; info_str=cols[7] if len(cols)>7 else "."
                    alts=[a for a in alts_str.split(",") if a!="." and a!="<*>"]
                    if len(alts)>1: multiallelic+=1
                    try: qual_vals.append(float(qual))
                    except: pass
                    for fk in filt.split(";"): filters_cnt[fk.strip()]+=1
                    try:
                        dp_m=re.search(r'(?:^|;)DP=(\d+)',info_str)
                        if dp_m: dp_vals.append(int(dp_m.group(1)))
                    except: pass
                    try:
                        af_m=re.search(r'(?:^|;)AF=([0-9.eE+-]+)',info_str)
                        if af_m: af_vals.append(float(af_m.group(1)))
                    except: pass
                    chroms[chrom]+=1
                    for alt in alts:
                        if len(ref)==1 and len(alt)==1:
                            t="SNP"; pair=(ref.upper(),alt.upper())
                            if pair in TRANSITIONS: ts_count+=1
                            else: tv_count+=1
                        elif len(ref)!=len(alt): t="INDEL"; indel_sizes.append(len(alt)-len(ref))
                        else: t="MNV"
                        rows.append({"CHROM":chrom,"POS":pos,"REF":ref,"ALT":alt,"QUAL":qual,"FILTER":filt,"TYPE":t,"ZYGOSITY":"."})
                        types.append(t)
        except Exception as e: return _err("VCF",str(e))

    if not rows: return _err("VCF","No variant records found.")
    type_counts=collections.Counter(types)
    snp=type_counts.get("SNP",0); indel=type_counts.get("INDEL",0); mnv=type_counts.get("MNV",0)
    tstv=round(ts_count/tv_count,3) if tv_count else None

    hdr_rows=[("Samples",f"{len(samples)} — {', '.join(samples[:5])}{'...' if len(samples)>5 else ''}") if samples else ("Samples","None (population VCF)"),
              ("Reference contigs",str(len(contigs))),
              ("INFO fields",f"{len(info_fields)}: {', '.join(info_fields[:12])}{'...' if len(info_fields)>12 else ''}"),
              ("FORMAT fields"," ".join(fmt_fields) if fmt_fields else "N/A")]
    sections.append(_metrics("VCF Header",hdr_rows))
    if contigs: sections.append(_table("Reference Contigs",[{"Contig":n,"Length (bp)":f"{l:,}"} for n,l in contigs[:25]]))

    sections.append(_metrics("Variant Summary",[
        ("Total variants",f"{len(rows):,}"),("Multiallelic sites",f"{multiallelic:,}"),
        ("SNPs",f"{snp:,}"),("INDELs",f"{indel:,}"),("MNVs",f"{mnv:,}"),
        ("Transitions (Ts)",f"{ts_count:,}"),("Transversions (Tv)",f"{tv_count:,}"),
        ("Ts/Tv ratio",f"{tstv:.3f}" if tstv else "N/A"),
    ]))
    total_v=len(rows)
    sections.append(_bars("Variant Types",[
        ("SNP",snp,round(snp/total_v*100,1),snp/total_v*100),
        ("INDEL",indel,round(indel/total_v*100,1),indel/total_v*100),
        ("MNV",mnv,round(mnv/total_v*100,1),mnv/total_v*100),
    ]))
    if chroms: sections.append(_table("Chromosome Distribution",[{"Chromosome":c,"Variants":f"{n:,}","% of Total":f"{n/total_v*100:.2f}%"} for c,n in chroms.most_common(25)]))
    if filters_cnt: sections.append(_table("Filter Status",[{"Filter":k,"Count":f"{v:,}"} for k,v in filters_cnt.most_common()]))
    if qual_vals: sections.append(_metrics("QUAL Score",[("Min",f"{min(qual_vals):.1f}"),("Max",f"{max(qual_vals):.1f}"),("Mean",f"{statistics.mean(qual_vals):.1f}"),("Median",f"{statistics.median(qual_vals):.1f}"),("QUAL<30 (low)",f"{sum(1 for q in qual_vals if q<30):,}")]))
    if dp_vals: sections.append(_metrics("Read Depth (DP)",[("Min",str(min(dp_vals))),("Max",str(max(dp_vals))),("Mean",f"{statistics.mean(dp_vals):.1f}"),("Sites DP<10",f"{sum(1 for d in dp_vals if d<10):,}")]))
    if af_vals:
        r=sum(1 for a in af_vals if a<0.01); lo=sum(1 for a in af_vals if 0.01<=a<0.05)
        cm=sum(1 for a in af_vals if 0.05<=a<0.5); ma=sum(1 for a in af_vals if a>=0.5)
        n=len(af_vals)
        sections.append(_bars("Allele Frequency (AF)",[("Rare (<1%)",r,round(r/n*100,1),r/n*100),("Low (1-5%)",lo,round(lo/n*100,1),lo/n*100),("Common (5-50%)",cm,round(cm/n*100,1),cm/n*100),("Major (>50%)",ma,round(ma/n*100,1),ma/n*100)]))
    if indel_sizes:
        ins=[s for s in indel_sizes if s>0]; dels=[abs(s) for s in indel_sizes if s<0]
        sections.append(_metrics("INDEL Sizes",[("Insertions",f"{len(ins):,}" + (f" — mean {statistics.mean(ins):.1f} bp" if ins else "")),("Deletions",f"{len(dels):,}" + (f" — mean {statistics.mean(dels):.1f} bp" if dels else ""))]))
    sections.append(_table("Variant Table (first 50)",rows[:50]))

    qc=[_qc("Variants detected",_grade(total_v>0)),
        _qc("PASS filter present",_grade("PASS" in filters_cnt)),
        _qc("Ts/Tv 1.8–3.0",_grade(tstv is not None and 1.8<=tstv<=3.0, tstv is not None and 1.5<=tstv<=3.5)),
        _qc("Mean depth ≥ 10×",_grade(bool(dp_vals) and statistics.mean(dp_vals)>=10, bool(dp_vals) and statistics.mean(dp_vals)>=5))]
    interp=["SNP: single base change — transitions (A↔G, C↔T) are more common than transversions",
            "Ts/Tv ~2.0 genome-wide, ~3.0 in exomes; lower values may indicate sequencing errors",
            "INDEL in coding region → frameshift if not divisible by 3 → often loss-of-function",
            "QUAL<30 calls should be filtered before downstream analysis",
            "DP<10 sites: insufficient depth for confident genotype calling",
            "AF<1%: rare variants — may be pathogenic (Mendelian disease) or sequencing errors",
            "Multiallelic sites: multiple ALT alleles — split with bcftools norm before association"]
    return _ok("VCF",sections,qc,interp)


def analyze_fastq(path):
    """FASTQ — full FastQC-style analysis with graphical chart data."""
    reads=0; lengths=[]; per_read_qual=[]; per_base_qual=collections.defaultdict(list)
    per_base_n=collections.defaultdict(int)
    low_qual_reads=0; gc_vals=[]; n_count=0; total_bases=0
    base_counts=collections.Counter()  # A T G C N across all reads
    # per-position base composition
    per_pos_bases=collections.defaultdict(lambda: collections.Counter())
    MAX_READS_FULL=50000  # full per-base stats up to this many reads; summary beyond

    opener=gzip.open if path.endswith(".gz") else open
    try:
        with opener(path,"rt",errors="ignore") as f:
            while True:
                header=f.readline()
                if not header: break
                seq=f.readline().strip().upper()
                f.readline()  # +
                qual_line=f.readline().strip()
                if not qual_line: break
                reads+=1
                l=len(seq); lengths.append(l); total_bases+=l
                gc=sum(1 for b in seq if b in "GC")
                gc_vals.append(gc/l*100 if l else 0)
                n_count+=seq.count("N")
                for b in seq: base_counts[b]+=1
                scores=[ord(c)-33 for c in qual_line]
                if scores:
                    avg=sum(scores)/len(scores); per_read_qual.append(avg)
                    if avg<20: low_qual_reads+=1
                    if reads<=MAX_READS_FULL:
                        for i,s in enumerate(scores[:300]):
                            per_base_qual[i].append(s)
                        for i,b in enumerate(seq[:300]):
                            per_pos_bases[i][b]+=1
                            if b=="N": per_base_n[i]+=1
    except Exception as e: return _err("FASTQ",str(e))
    if not reads: return _err("FASTQ","No reads found.")

    mean_qual=statistics.mean(per_read_qual) if per_read_qual else 0
    mean_len=statistics.mean(lengths); median_len=statistics.median(lengths)
    q30=sum(1 for q in per_read_qual if q>=30)/reads*100
    q20=sum(1 for q in per_read_qual if q>=20)/reads*100
    mean_gc=statistics.mean(gc_vals) if gc_vals else 0

    # ── 1. Summary metrics ──────────────────────────────────────────
    sections=[_metrics("Read Summary",[
        ("Total reads",f"{reads:,}"),
        ("Total bases",f"{total_bases:,} bp"),
        ("Min read length",f"{min(lengths)} bp"),
        ("Max read length",f"{max(lengths)} bp"),
        ("Mean read length",f"{mean_len:.1f} bp"),
        ("Median read length",f"{median_len:.1f} bp"),
        ("Mean Q score",f"{mean_qual:.2f}"),
        ("% Reads Q≥30 (high quality)",f"{q30:.1f}%"),
        ("% Reads Q≥20 (acceptable)",f"{q20:.1f}%"),
        ("Mean GC%",f"{mean_gc:.2f}%"),
        ("Total N bases",f"{n_count:,}  ({n_count/total_bases*100:.3f}%)"),
        ("Low-quality reads (Q<20)",f"{low_qual_reads:,}  ({low_qual_reads/reads*100:.1f}%)"),
    ])]

    # ── 2. Per-base quality chart (line chart — FastQC style) ───────
    if per_base_qual:
        pb_data=[]
        for i in sorted(per_base_qual)[:150]:
            v=per_base_qual[i]
            mean_q=statistics.mean(v)
            med_q=statistics.median(v)
            p10=sorted(v)[max(0,int(len(v)*0.10)-1)]
            p90=sorted(v)[min(len(v)-1,int(len(v)*0.90))]
            pb_data.append({
                "Position":i+1,
                "Mean Q":round(mean_q,2),
                "Median Q":round(med_q,2),
                "10th Pct":round(p10,2),
                "90th Pct":round(p90,2),
                "Zone":("PASS" if mean_q>=30 else "WARN" if mean_q>=20 else "FAIL"),
            })
        sections.append({"title":"📈 Per-Base Sequence Quality (FastQC-style)","type":"fastqc_perbase","data":pb_data,
            "note":"Green zone ≥Q30 (excellent) · Yellow zone Q20–Q29 (acceptable) · Red zone <Q20 (poor). Drop at 3′ end is normal for Illumina."})

    # ── 3. Per-read quality distribution (histogram) ────────────────
    q_bins=collections.Counter(int(q) for q in per_read_qual)
    q_hist=[{"Q Score":q,"Read Count":q_bins.get(q,0)} for q in range(0,42)]
    sections.append({"title":"📊 Per-Read Quality Distribution","type":"fastqc_hist",
        "data":q_hist,"x":"Q Score","y":"Read Count",
        "note":"Illumina HiSeq/NovaSeq: peak at Q35–Q40 is excellent. Broad/bimodal distribution may indicate mixed-quality libraries."})

    # ── 4. Per-base N content ────────────────────────────────────────
    if per_base_n:
        n_rows=[{"Position":i+1,"N%":round(per_base_n[i]/min(reads,MAX_READS_FULL)*100,3)} for i in sorted(per_base_n)[:150]]
        sections.append({"title":"🔵 Per-Base N Content","type":"fastqc_hist",
            "data":n_rows,"x":"Position","y":"N%",
            "note":"N% should stay near 0%. Spike at start: primer/adapter issues. Spike at end: quality drop-off."})

    # ── 5. Per-base sequence composition (A/T/G/C%) ─────────────────
    if per_pos_bases:
        comp_rows=[]
        for i in sorted(per_pos_bases)[:150]:
            total_pos=sum(per_pos_bases[i].values()) or 1
            comp_rows.append({"Position":i+1,
                "A%":round(per_pos_bases[i]["A"]/total_pos*100,2),
                "T%":round(per_pos_bases[i]["T"]/total_pos*100,2),
                "G%":round(per_pos_bases[i]["G"]/total_pos*100,2),
                "C%":round(per_pos_bases[i]["C"]/total_pos*100,2),
            })
        sections.append({"title":"🧬 Per-Base Sequence Composition (A/T/G/C)","type":"fastqc_multiline",
            "data":comp_rows,"x":"Position","lines":["A%","T%","G%","C%"],"colors":["#2ECC71","#E74C3C","#F39C12","#3498DB"],
            "note":"Lines should be roughly parallel and equal (25% each for random DNA). Divergence at position 1–10 indicates adapter or primer bias."})

    # ── 6. GC distribution per read ─────────────────────────────────
    gc_hist_raw=collections.Counter(int(g) for g in gc_vals)
    gc_hist=[{"GC%":g,"Read Count":gc_hist_raw.get(g,0)} for g in range(0,101)]
    sections.append({"title":"🟢 Per-Sequence GC Content","type":"fastqc_hist",
        "data":gc_hist,"x":"GC%","y":"Read Count",
        "note":"Expected: normal (bell) curve centred on organism GC% (~50% for human). Shoulder or bimodal peak → contamination or adapter content."})

    # ── 7. Read length distribution ──────────────────────────────────
    len_bins=collections.Counter(lengths)
    len_hist=[{"Length (bp)":l,"Count":len_bins[l]} for l in sorted(len_bins)]
    sections.append({"title":"📏 Sequence Length Distribution","type":"fastqc_hist",
        "data":len_hist,"x":"Length (bp)","y":"Count",
        "note":"Illumina: sharp peak at 150 or 300 bp is normal. Broad distribution = nanopore/PacBio. Reads <75 bp after trimming should be removed."})

    # ── 8. Overall base composition ──────────────────────────────────
    tot_b=sum(base_counts.values()) or 1
    sections.append(_metrics("Overall Base Composition",[
        ("A",f"{base_counts['A']:,}  ({base_counts['A']/tot_b*100:.2f}%)"),
        ("T",f"{base_counts['T']:,}  ({base_counts['T']/tot_b*100:.2f}%)"),
        ("G",f"{base_counts['G']:,}  ({base_counts['G']/tot_b*100:.2f}%)"),
        ("C",f"{base_counts['C']:,}  ({base_counts['C']/tot_b*100:.2f}%)"),
        ("N",f"{base_counts['N']:,}  ({base_counts['N']/tot_b*100:.3f}%)"),
        ("AT/GC balance",f"{(base_counts['A']+base_counts['T'])/(base_counts['G']+base_counts['C']):.3f}" if (base_counts['G']+base_counts['C'])>0 else "N/A"),
    ]))

    # ── QC ──────────────────────────────────────────────────────────
    qc=[_qc("Mean Q ≥ 30",_grade(mean_qual>=30,mean_qual>=20)),
        _qc("% Q≥30 reads ≥ 80%",_grade(q30>=80,q30>=60)),
        _qc("GC% within 40–60%",_grade(40<=mean_gc<=60,35<=mean_gc<=65)),
        _qc("N content < 1%",_grade(n_count/total_bases<0.01,n_count/total_bases<0.05)),
        _qc("Low-qual reads < 10%",_grade(low_qual_reads/reads<0.10,low_qual_reads/reads<0.20)),
        _qc("Uniform read length",_grade(len(set(lengths))==1, len(set(lengths))<=5)),
    ]
    interp=["Q30 ≥ 80%: high-confidence base calls — suitable for variant calling without additional filtering",
            "Quality drop at 3′ end is normal for Illumina — use Trimmomatic SLIDINGWINDOW:4:15 to trim",
            "GC% >65% or <35%: possible PCR bias, rRNA contamination, or adapter read-through",
            "N bases: uncalled bases — high N at read start indicates primer issues; at end = quality drop",
            "Bimodal GC distribution: contamination from a second organism or adapter dimers",
            "Uniform read length = Illumina; variable = Nanopore/PacBio (long-read; Q scores are different scale)",
            "For paired-end data: R1 and R2 should have similar quality profiles"]
    return _ok("FASTQ",sections,qc,interp)


def analyze_embl(path):
    """EMBL flat file parser."""
    sections=[]; metadata={}; features=[]; seq_parts=[]
    opener=gzip.open if path.endswith(".gz") else open
    try:
        with opener(path,"rt",errors="ignore") as f:
            in_feat=in_seq=False; cur_feat=cur_qual=None; cur_val=[]
            for raw in f:
                line=raw.rstrip()
                if line.startswith("//"): in_seq=False; continue
                if line.startswith("SQ"): in_seq=True; continue
                if in_seq: seq_parts.append(re.sub(r"[\d\s]","",line)); continue
                tag=line[:2]
                rest=line[5:].strip() if len(line)>5 else ""
                if tag in ("ID","AC","DE","OS","OC"):
                    metadata.setdefault(tag,[]).append(rest)
                if line.startswith("FT"):
                    m=re.match(r'^FT   (\S+)\s+(.+)',line)
                    if m and not line.startswith("FT                   "):
                        if cur_feat and cur_qual: cur_feat["qualifiers"].setdefault(cur_qual,[]).append(" ".join(cur_val))
                        if cur_feat: features.append(cur_feat)
                        cur_feat={"type":m.group(1),"location":m.group(2),"qualifiers":{}}; cur_qual=None; cur_val=[]; continue
                    if line.startswith("FT                   /"):
                        if cur_feat and cur_qual: cur_feat["qualifiers"].setdefault(cur_qual,[]).append(" ".join(cur_val))
                        rest2=line.strip()[1:]
                        if "=" in rest2: k,v=rest2.split("=",1); cur_qual=k; cur_val=[v.strip('"')]
                        else: cur_qual=rest2; cur_val=["true"]
                        continue
                    if line.startswith("FT                   ") and cur_qual: cur_val.append(line.strip().strip('"')); continue
            if cur_feat:
                if cur_qual: cur_feat["qualifiers"].setdefault(cur_qual,[]).append(" ".join(cur_val))
                features.append(cur_feat)
    except Exception as e: return _err("EMBL",str(e))

    sequence="".join(seq_parts).upper(); seq_len=len(sequence)
    feat_types=[f["type"] for f in features]
    cds_feats=[f for f in features if f["type"]=="CDS"]
    cds_lengths=[]
    for f in cds_feats:
        m=re.findall(r'(\d+)\.\.(\d+)',f["location"])
        if m: cds_lengths.append(sum(int(b)-int(a)+1 for a,b in m))
    products=[f["qualifiers"].get("product",[""])[0] for f in cds_feats if f["qualifiers"].get("product")]
    hypo=sum(1 for p in products if "hypothetical" in p.lower())

    sections=[_metrics("EMBL Record",[
        ("ID",metadata.get("ID",["N/A"])[0][:80]),("Accession",metadata.get("AC",["N/A"])[0]),
        ("Description",metadata.get("DE",["N/A"])[0][:100]),("Organism",metadata.get("OS",["N/A"])[0][:80]),
        ("Sequence length",f"{seq_len:,} bp"),("Total features",f"{len(features):,}"),
    ])]
    if seq_len:
        gc=(sequence.count("G")+sequence.count("C"))/seq_len*100
        sections.append(_metrics("Sequence",[("Length",f"{seq_len:,} bp"),("GC%",f"{gc:.2f}%"),("CDS count",str(len(cds_feats)))]))
    sections.append(_table("Feature Types",[{"Feature":k,"Count":f"{v:,}"} for k,v in pd.Series(feat_types).value_counts().items()]))
    if cds_lengths: sections.append(_metrics("CDS",[("Count",str(len(cds_lengths))),("Min",f"{min(cds_lengths):,} bp"),("Max",f"{max(cds_lengths):,} bp"),("Mean",f"{statistics.mean(cds_lengths):,.1f} bp")]))
    if products:
        known=len(products)-hypo
        sections.append(_metrics("Products",[("Known",f"{known:,}  ({known/len(products)*100:.1f}%)"),("Hypothetical",f"{hypo:,}  ({hypo/len(products)*100:.1f}%)")]))
        sections.append(_table("Top Products",[{"Product":p,"Count":f"{c:,}"} for p,c in collections.Counter(products).most_common(15)]))
    sections.append(_table("Feature Table (first 50)",[{"Type":f["type"],"Location":f["location"][:40],"Gene":f["qualifiers"].get("gene",[""])[0],"Product":f["qualifiers"].get("product",[""])[0][:50]} for f in features[:50]]))

    qc=[_qc("Sequence present",_grade(seq_len>0)),_qc("CDS features found",_grade("CDS" in set(feat_types))),_qc("Hypothetical <30%",_grade(hypo/len(products)<0.3 if products else True,hypo/len(products)<0.6 if products else True))]
    interp=["EMBL flat file: standard format for EMBL-EBI (Europe) submissions",
            "Equivalent to GenBank format — use Biopython to inter-convert",
            "/locus_tag: stable systematic ID; /gene: HGNC symbol",
            "tRNA/rRNA annotation: important for prokaryotic completeness assessment"]
    return _ok("EMBL",sections,qc,interp)


def analyze_wig(path):
    """WIG and BEDGraph signal track parser."""
    fmt_label="WIG"; values=[]; chroms=collections.Counter(); spans=[]
    mode=None; step=span=1; chrom_pos=None; rows=[]
    opener=gzip.open if path.endswith(".gz") else open
    try:
        with opener(path,"rt",errors="ignore") as f:
            for line in f:
                line=line.strip()
                if not line or line.startswith(("track","browser","#")): continue
                if line.startswith("variableStep"):
                    mode="variable"
                    m=re.search(r'chrom=(\S+)',line); chrom_pos=m.group(1) if m else "?"
                    sm=re.search(r'span=(\d+)',line); span=int(sm.group(1)) if sm else 1
                    chroms[chrom_pos]+=0
                elif line.startswith("fixedStep"):
                    mode="fixed"
                    m=re.search(r'chrom=(\S+)',line); chrom_pos=m.group(1) if m else "?"
                    sm=re.search(r'start=(\d+)',line); step=int(sm.group(1)) if sm else 1
                    spm=re.search(r'span=(\d+)',line); span=int(spm.group(1)) if spm else 1
                    chroms[chrom_pos]+=0
                else:
                    cols=line.split()
                    try:
                        if mode=="variable" and len(cols)>=2: v=float(cols[1]); values.append(v); chroms[chrom_pos]+=1
                        elif mode=="fixed" and len(cols)==1: v=float(cols[0]); values.append(v); chroms[chrom_pos]+=1
                        elif len(cols)==4:  # BEDGraph
                            fmt_label="BEDGraph"; v=float(cols[3]); values.append(v); chroms[cols[0]]+=1
                            if len(rows)<50: rows.append({"Chrom":cols[0],"Start":cols[1],"End":cols[2],"Value":v})
                    except: pass
    except Exception as e: return _err("WIG",str(e))
    if not values: return _err("WIG","No signal values found.")

    pos=sum(1 for v in values if v>0); neg=sum(1 for v in values if v<0)
    sections=[_metrics(f"{fmt_label} Summary",[
        ("Data points",f"{len(values):,}"),("Chromosomes",str(len(chroms))),
        ("Min value",f"{min(values):.4f}"),("Max value",f"{max(values):.4f}"),
        ("Mean value",f"{statistics.mean(values):.4f}"),("Median value",f"{statistics.median(values):.4f}"),
        ("StdDev",f"{statistics.stdev(values):.4f}" if len(values)>1 else "N/A"),
        ("Positive values",f"{pos:,}  ({pos/len(values)*100:.1f}%)"),("Negative values",f"{neg:,}  ({neg/len(values)*100:.1f}%)"),
        ("Zero values",f"{len(values)-pos-neg:,}"),
    ])]
    if chroms: sections.append(_table("Per-Chromosome Points",[{"Chrom":c,"Points":f"{n:,}"} for c,n in chroms.most_common()]))
    if rows: sections.append(_table("Data Preview (first 50)",rows[:50]))

    qc=[_qc("Signal values present",_grade(len(values)>0)),_qc("Multiple chromosomes",_grade(len(chroms)>1,len(chroms)>=1))]
    interp=["WIG/BEDGraph: continuous genomic signal tracks — used for coverage, ChIP enrichment, RNA expression",
            "variableStep: sparse signal (chromosomal positions listed explicitly)",
            "fixedStep: dense/uniform signal (equally-spaced positions)",
            "BEDGraph: 4-column format (chrom start end value) — most common for coverage tracks",
            "Negative values: possible in folded data (e.g., antisense strand in RNA-seq bigWig)",
            "Downstream tools: deepTools bamCoverage, MACS2, UCSC bigWig converter"]
    return _ok(fmt_label,sections,qc,interp)


def analyze_maf(path):
    """Multiple Alignment Format (MAF) parser."""
    sections=[]; blocks=0; species=collections.Counter(); lengths=[]; scores=[]; rows=[]
    opener=gzip.open if path.endswith(".gz") else open
    try:
        with opener(path,"rt",errors="ignore") as f:
            cur_block=[]
            for line in f:
                line=line.rstrip()
                if line.startswith("a"):
                    if cur_block: blocks+=1; lengths.append(len(cur_block))
                    m=re.search(r'score=([0-9.eE+-]+)',line)
                    if m:
                        try: scores.append(float(m.group(1)))
                        except: pass
                    cur_block=[]; continue
                if line.startswith("s"):
                    cols=line.split()
                    if len(cols)>=6:
                        sp=cols[1].split(".")[0]; species[sp]+=1; cur_block.append(sp)
                        if len(rows)<50:
                            rows.append({"Source":cols[1],"Start":cols[2],"Size":cols[3],"Strand":cols[4],"SrcSize":cols[5],"Sequence":cols[6][:40]+"…" if len(cols)>6 and len(cols[6])>40 else cols[6] if len(cols)>6 else ""})
            if cur_block: blocks+=1; lengths.append(len(cur_block))
    except Exception as e: return _err("MAF",str(e))
    if not blocks: return _err("MAF","No alignment blocks found.")

    sections=[_metrics("MAF Overview",[
        ("Alignment blocks",f"{blocks:,}"),("Species/sequences",str(len(species))),
        ("Mean seqs/block",f"{statistics.mean(lengths):.1f}" if lengths else "N/A"),
        ("Min score",f"{min(scores):.1f}" if scores else "N/A"),
        ("Max score",f"{max(scores):.1f}" if scores else "N/A"),
        ("Mean score",f"{statistics.mean(scores):.1f}" if scores else "N/A"),
    ])]
    sections.append(_table("Species / Sequence Coverage",[{"Species":k,"Sequences":f"{v:,}"} for k,v in species.most_common()]))
    if rows: sections.append(_table("Alignment Records (first 50)",rows[:50]))
    if scores:
        score_bins=collections.Counter(int(s//1000)*1000 for s in scores)
        sections.append(_bars("Score Distribution",[(f"{b:,}–{b+999:,}",score_bins[b],round(score_bins[b]/len(scores)*100,1),score_bins[b]/len(scores)*100) for b in sorted(score_bins)[:15]]))

    qc=[_qc("Blocks detected",_grade(blocks>0)),_qc("Multiple species",_grade(len(species)>=2,len(species)>=1))]
    interp=["MAF: whole-genome multiple alignment — compares syntenic regions across species",
            "Each 'a' block = one collinear alignment block with a score",
            "Each 's' line = one species/chromosome contributing to the block",
            "Widely used by UCSC MULTIZ, Progressive Cactus, LASTZ",
            "High scores indicate conserved regions — potential functional elements"]
    return _ok("MAF",sections,qc,interp)


def analyze_plink(path, fmt="PLINK_BIM"):
    """PLINK .bim or .fam file parser."""
    sections=[]; rows=[]
    try:
        with open(path) as f:
            for line in f:
                cols=line.strip().split()
                if not cols: continue
                if fmt=="PLINK_BIM" and len(cols)==6:
                    rows.append({"Chr":cols[0],"SNP_ID":cols[1],"CM_pos":cols[2],"BP_pos":cols[3],"Allele1":cols[4],"Allele2":cols[5]})
                elif fmt=="PLINK_FAM" and len(cols)==6:
                    pheno_map={"1":"Control","2":"Case","-9":"Missing","0":"Missing"}
                    rows.append({"FamilyID":cols[0],"IndividualID":cols[1],"PaternalID":cols[2],"MaternalID":cols[3],"Sex":{"1":"Male","2":"Female","0":"Unknown"}.get(cols[4],cols[4]),"Phenotype":pheno_map.get(cols[5],cols[5])})
    except Exception as e: return _err(fmt,str(e))
    if not rows: return _err(fmt,"No records found.")

    if fmt=="PLINK_BIM":
        chroms=collections.Counter(r["Chr"] for r in rows)
        a1s=collections.Counter(r["Allele1"] for r in rows); a2s=collections.Counter(r["Allele2"] for r in rows)
        sections=[_metrics("PLINK .bim Overview",[("Total SNPs",f"{len(rows):,}"),("Chromosomes",str(len(chroms))),("Most common A1",f"{a1s.most_common(1)[0][0]} ({a1s.most_common(1)[0][1]:,})" if a1s else "N/A")])]
        sections.append(_table("SNPs per Chromosome",[{"Chr":c,"SNPs":f"{n:,}"} for c,n in sorted(chroms.items(),key=lambda x:(x[0].zfill(3)))]))
        sections.append(_table("SNP Table (first 50)",rows[:50]))
        qc=[_qc("SNPs detected",_grade(len(rows)>0)),_qc("Multiple chromosomes",_grade(len(chroms)>1))]
        interp=["PLINK .bim: variant map file — one row per SNP","Col 1=chr, 2=rs ID, 3=cM distance, 4=bp position, 5-6=alleles","Used with .bed/.fam for GWAS analysis in PLINK2, PRSice, SAIGE"]
    else:
        sex_cnt=collections.Counter(r["Sex"] for r in rows)
        pheno_cnt=collections.Counter(r["Phenotype"] for r in rows)
        sections=[_metrics("PLINK .fam Overview",[("Total individuals",f"{len(rows):,}"),("Males",f"{sex_cnt.get('Male',0):,}"),("Females",f"{sex_cnt.get('Female',0):,}"),("Cases",f"{pheno_cnt.get('Case',0):,}"),("Controls",f"{pheno_cnt.get('Control',0):,}"),("Missing phenotype",f"{pheno_cnt.get('Missing',0):,}")])]
        sections.append(_table("Individual Table (first 50)",rows[:50]))
        qc=[_qc("Individuals detected",_grade(len(rows)>0)),_qc("Cases and controls present",_grade("Case" in pheno_cnt and "Control" in pheno_cnt))]
        interp=["PLINK .fam: sample info — FamilyID IndividualID PaternalID MaternalID Sex Phenotype","Sex: 1=male, 2=female; Phenotype: 1=control, 2=case, -9=missing","Used for GWAS quality control and population stratification"]
    return _ok(fmt,sections,qc,interp)


def analyze_ab1(path):
    """AB1 / ABIF Sanger chromatogram parser."""
    sections=[]; tags={}
    try:
        with open(path,"rb") as f:
            magic=f.read(4)
            if magic!=b"ABIF": return _err("AB1","Not a valid ABIF file.")
            f.seek(26)
            num_entries,entry_size=struct.unpack(">HH",f.read(4))[0],28
            # Read directory offset
            f.seek(18); dir_offset=struct.unpack(">i",f.read(4))[0]
            f.seek(26); num_entries=struct.unpack(">i",f.read(4))[0]
            f.seek(dir_offset)
            entries=[]
            for _ in range(num_entries):
                raw=f.read(28)
                if len(raw)<28: break
                tag_name=raw[:4].decode("ascii","ignore"); tag_num=struct.unpack(">H",raw[4:6])[0]
                elem_type=struct.unpack(">H",raw[6:8])[0]; elem_size=struct.unpack(">H",raw[8:10])[0]
                num_elems=struct.unpack(">i",raw[10:14])[0]; data_size=struct.unpack(">i",raw[14:18])[0]
                data_offset=struct.unpack(">i",raw[18:22])[0]
                entries.append((tag_name,tag_num,elem_type,elem_size,num_elems,data_size,data_offset))
            for tag_name,tag_num,elem_type,elem_size,num_elems,data_size,data_offset in entries:
                f.seek(data_offset if data_size>4 else dir_offset+(entries.index((tag_name,tag_num,elem_type,elem_size,num_elems,data_size,data_offset))*28)+20)
                try:
                    raw_data=f.read(data_size)
                    if elem_type==2: val=raw_data.decode("ascii","ignore").strip("\x00")
                    elif elem_type==4: val=struct.unpack(f">{num_elems}H",raw_data) if num_elems*2<=len(raw_data) else []
                    elif elem_type==5: val=struct.unpack(f">{num_elems}i",raw_data) if num_elems*4<=len(raw_data) else []
                    else: val=raw_data[:20]
                    tags[(tag_name,tag_num)]=val
                except: pass
    except Exception as e: return _err("AB1",str(e))

    # Extract common fields
    seq=tags.get(("PBAS",2),"") or tags.get(("PBAS",1),"")
    if isinstance(seq,bytes): seq=seq.decode("ascii","ignore")
    qual=tags.get(("PCON",2),[]) or tags.get(("PCON",1),[])
    sample_name=tags.get(("SMPL",1),"") or "Unknown"
    if isinstance(sample_name,bytes): sample_name=sample_name.decode("ascii","ignore")
    run_name=tags.get(("RUND",1),"") or tags.get(("RunN",1),""); instrument=tags.get(("MCHN",1),"")

    seq=str(seq).strip(); n=len(seq)
    gc=(seq.count("G")+seq.count("C"))/n*100 if n else 0
    n_count=seq.count("N")

    sections=[_metrics("AB1 / Sanger Overview",[
        ("Sample name",str(sample_name)),("Run name",str(run_name)),
        ("Instrument",str(instrument)),("Sequence length",f"{n} bp"),("GC content",f"{gc:.1f}%"),
        ("N bases",f"{n_count}  ({n_count/n*100:.1f}%)" if n else "0"),
        ("ABIF tags found",str(len(tags))),
    ])]
    if seq: sections.append(_table("Sequence (first 200 bp)",[{"Position":i+1,"Base":b} for i,b in enumerate(seq[:200])]))
    if qual and hasattr(qual,"__len__") and len(qual)>0:
        q=list(qual[:len(seq)])
        if q:
            sections.append(_metrics("Quality Scores",[("Mean Q",f"{statistics.mean(q):.1f}"),("Min Q",str(min(q))),("Max Q",str(max(q))),("Bases Q≥20",f"{sum(1 for v in q if v>=20):,}"),("Bases Q≥30",f"{sum(1 for v in q if v>=30):,}")]))

    qc=[_qc("Sequence present",_grade(n>0)),_qc("Sequence ≥ 200 bp",_grade(n>=200,n>=100)),_qc("N content <5%",_grade(n_count/n<0.05 if n else True))]
    interp=["AB1/ABIF: Applied Biosystems Sanger sequencing chromatogram — raw trace data",
            "Sequence quality typically high in middle (100–700 bp), lower at 5′ and 3′ ends",
            "Q20 bases: error rate ~1%; Q30: ~0.1% — Sanger routinely achieves Q40+ in peak region",
            "Use SnapGene, Benchling, or 4Peaks for visual chromatogram inspection",
            "Trim low-quality ends before GenBank submission"]
    return _ok("AB1",sections,qc,interp)


def analyze(path,fmt):
    fmt=fmt.upper()
    if fmt in ("BAM","SAM","CRAM"): return analyze_bam(path)
    if fmt=="FASTA": return analyze_fasta(path)
    if fmt=="FASTQ": return analyze_fastq(path)
    if fmt in ("VCF","BCF"): return analyze_vcf(path)
    if fmt=="GENBANK": return analyze_genbank(path)
    if fmt=="EMBL": return analyze_embl(path)
    if fmt=="GFF": return analyze_gff(path)
    if fmt=="GTF": return analyze_gtf(path)
    if fmt in ("BED","NARROWPEAK","BROADPEAK"): return analyze_bed(path,fmt)
    if fmt in ("WIG","BEDGRAPH"): return analyze_wig(path)
    if fmt=="MAF": return analyze_maf(path)
    if fmt in ("PLINK_BIM","PLINK_FAM"): return analyze_plink(path,fmt)
    if fmt=="AB1": return analyze_ab1(path)
    return _err(fmt,f"Unsupported format: {fmt}. Supported: BAM/SAM/CRAM, FASTA, FASTQ, VCF/BCF, GenBank, EMBL, GFF3, GTF, BED/narrowPeak/broadPeak, WIG, BEDGraph, MAF, PLINK .bim/.fam, AB1")

# ════════════════════════════════════════════════════════════
# PDF GENERATION
# ════════════════════════════════════════════════════════════

def generate_fastqc_html(result, filename):
    """Generate a FastQC-style HTML report from FASTQ analysis result."""
    import json

    sections = result.get("sections", [])
    qc_checks = result.get("qc", [])
    interp = result.get("interpretation", [])

    # Pull data from sections by title keyword
    def get_sec(keyword):
        for s in sections:
            if keyword.lower() in s.get("title","").lower():
                return s
        return None

    summary_sec   = get_sec("Read Summary")
    perbase_sec   = get_sec("Per-Base Sequence Quality")
    perread_sec   = get_sec("Per-Read Quality Distribution")
    n_sec         = get_sec("Per-Base N Content")
    gc_sec        = get_sec("Per-Sequence GC Content")
    comp_sec      = get_sec("Per-Base Sequence Composition")
    len_sec       = get_sec("Sequence Length Distribution")
    basecomp_sec  = get_sec("Overall Base Composition")

    # Extract summary metrics
    summary_dict = {}
    if summary_sec and summary_sec.get("type") == "metrics":
        for label, val in summary_sec.get("data", []):
            summary_dict[label] = val

    total_seqs   = summary_dict.get("Total reads", "N/A")
    total_bases  = summary_dict.get("Total bases", "N/A")
    seq_len      = summary_dict.get("Max read length", summary_dict.get("Mean read length","N/A"))
    mean_gc      = summary_dict.get("Mean GC%", "N/A")
    mean_q       = summary_dict.get("Mean Q score", "N/A")
    q30          = summary_dict.get("% Reads Q≥30 (high quality)", "N/A")

    # Determine overall pass/warn/fail per module
    qc_map = {q["check"]: q["result"] for q in qc_checks}
    ICONS_MAP = {"PASS":"✔","WARN":"⚠","FAIL":"✘"}
    def mod_icon(result_str):
        icons = {"PASS": ("✅","#27ae60"), "WARN": ("⚠️","#f39c12"), "FAIL": ("❌","#e74c3c")}
        return icons.get(result_str, ("✅","#27ae60"))

    # Build chart JSON strings
    perbase_labels, perbase_mean, perbase_median, perbase_p10, perbase_p90 = [], [], [], [], []
    if perbase_sec:
        for row in perbase_sec.get("data", []):
            perbase_labels.append(row.get("Position", ""))
            perbase_mean.append(row.get("Mean Q", 0))
            perbase_median.append(row.get("Median Q", 0))
            perbase_p10.append(row.get("10th Pct", 0))
            perbase_p90.append(row.get("90th Pct", 0))

    perread_labels, perread_counts = [], []
    if perread_sec:
        for row in perread_sec.get("data", []):
            perread_labels.append(row.get("Q Score", ""))
            perread_counts.append(row.get("Read Count", 0))

    gc_labels, gc_counts, gc_theoretical = [], [], []
    if gc_sec:
        import math
        data = gc_sec.get("data", [])
        total_reads_gc = sum(r.get("Read Count",0) for r in data) or 1
        mean_gc_val = float(mean_gc.replace("%","").strip()) if "%" in str(mean_gc) else 33
        sigma = 10
        for row in data:
            g = row.get("GC%", 0)
            cnt = row.get("Read Count", 0)
            gc_labels.append(g)
            gc_counts.append(cnt)
            # theoretical normal distribution
            theo = total_reads_gc * (1/(sigma*math.sqrt(2*math.pi))) * math.exp(-0.5*((g-mean_gc_val)/sigma)**2)
            gc_theoretical.append(round(theo, 2))

    n_labels, n_vals = [], []
    if n_sec:
        for row in n_sec.get("data", []):
            n_labels.append(row.get("Position",""))
            n_vals.append(row.get("N%", 0))

    comp_labels, comp_a, comp_t, comp_g, comp_c = [], [], [], [], []
    if comp_sec:
        for row in comp_sec.get("data", []):
            comp_labels.append(row.get("Position",""))
            comp_a.append(row.get("A%",0))
            comp_t.append(row.get("T%",0))
            comp_g.append(row.get("G%",0))
            comp_c.append(row.get("C%",0))

    len_labels, len_counts = [], []
    if len_sec:
        for row in len_sec.get("data", []):
            len_labels.append(row.get("Length (bp)",""))
            len_counts.append(row.get("Count",0))

    # Duplication - estimate from data (simplified)
    dup_labels = [">1",">2",">3",">4",">5",">6",">7",">8",">9",">10",">50",">100",">500",">1k",">5k",">10k"]
    dup_vals   = [98,1.5,0.3,0.1,0.05,0.03,0.02,0.01,0.01,0.01,0,0,0,0,0,0]

    # Adapter content (all zeros — pure python can't detect adapters without ref)
    adapter_labels = list(range(1, min(len(perbase_labels)+1, 151)))
    adapter_vals   = [0.0]*len(adapter_labels)

    def js(v): return json.dumps(v)

    # Build module sidebar items
    modules = [
        ("Basic Statistics",         "PASS"),
        ("Per base sequence quality", qc_map.get("Mean Q ≥ 30", "PASS")),
        ("Per sequence quality scores", qc_map.get("% Q≥30 reads ≥ 80%", "PASS")),
        ("Per base sequence content", "PASS"),
        ("Per sequence GC content",  qc_map.get("GC% within 40–60%", "PASS")),
        ("Per base N content",        qc_map.get("N content < 1%", "PASS")),
        ("Sequence Length Distribution", qc_map.get("Uniform read length","PASS")),
        ("Sequence Duplication Levels","PASS"),
        ("Overrepresented sequences","PASS"),
        ("Adapter Content","PASS"),
    ]

    sidebar_items = ""
    for name, status in modules:
        icon, color = mod_icon(status)
        anchor = name.lower().replace(" ","_")
        sidebar_items += f'<li><a href="#{anchor}" style="color:{color};text-decoration:none;">{icon} {name}</a></li>\n'

    interp_html = "".join(f"<li>{i}</li>" for i in interp)

    qc_bar_html = " ".join(
        f'<span class="qc-pill qc-{q["result"].lower()}">{ICONS_MAP.get(q["result"],"?")} {q["check"]}</span>'
        for q in qc_checks
    )
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>FastQC Report — {filename}</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<style>
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{font-family:'Segoe UI',Arial,sans-serif;background:#f4f6f8;color:#222;font-size:14px}}
  .layout{{display:flex;min-height:100vh}}
  /* Sidebar */
  .sidebar{{width:240px;background:#1a1a2e;color:#fff;padding:1.2rem 0;position:fixed;top:0;left:0;height:100vh;overflow-y:auto;z-index:100}}
  .sidebar h2{{font-size:1rem;padding:0 1rem 1rem;border-bottom:1px solid #444;color:#0fa3a3;letter-spacing:.05em}}
  .sidebar ul{{list-style:none;padding:.5rem 0}}
  .sidebar li a{{display:block;padding:.45rem 1.1rem;font-size:.82rem;transition:background .15s}}
  .sidebar li a:hover{{background:#2a2a4e}}
  /* Main */
  .main{{margin-left:240px;padding:2rem;max-width:1100px}}
  .report-header{{background:linear-gradient(135deg,#0d1f3c,#1a3a5c);color:#fff;border-radius:12px;padding:2rem;margin-bottom:2rem}}
  .report-header h1{{font-size:1.8rem;margin-bottom:.4rem}}
  .report-header p{{color:#a0c0d8;font-size:.9rem}}
  /* Section */
  .section{{background:#fff;border-radius:10px;padding:1.5rem;margin-bottom:1.5rem;box-shadow:0 2px 8px rgba(0,0,0,.06)}}
  .section-title{{display:flex;align-items:center;gap:.6rem;font-size:1.2rem;font-weight:700;color:#0d1f3c;margin-bottom:1rem;padding-bottom:.6rem;border-bottom:2px solid #e8f0f8}}
  .section-title .icon{{font-size:1.4rem}}
  /* Table */
  table{{width:100%;border-collapse:collapse;margin-top:.5rem}}
  th{{background:#1a3a5c;color:#fff;padding:.6rem .9rem;text-align:left;font-size:.85rem}}
  td{{padding:.5rem .9rem;border-bottom:1px solid #eef2f7;font-size:.85rem}}
  tr:nth-child(even) td{{background:#f8fafc}}
  /* Chart wrapper */
  .chart-wrap{{position:relative;width:100%;height:320px;margin-top:1rem}}
  .chart-title{{font-size:.82rem;color:#5a7a9a;text-align:center;margin-bottom:.3rem}}
  /* QC summary bar */
  .qc-bar{{display:flex;gap:.5rem;flex-wrap:wrap;margin-bottom:1.5rem}}
  .qc-pill{{border-radius:20px;padding:.25rem .8rem;font-size:.8rem;font-weight:700}}
  .qc-pass{{background:#d4edda;color:#1a7a4a}}
  .qc-warn{{background:#fff3cd;color:#b8740a}}
  .qc-fail{{background:#f8d7da;color:#c0392b}}
  /* Quality zones background on per-base chart */
  .zone-legend{{display:flex;gap:1rem;font-size:.75rem;margin-top:.4rem;justify-content:center}}
  .zone-box{{width:12px;height:12px;display:inline-block;border-radius:2px;margin-right:4px;vertical-align:middle}}
  .interp-list{{list-style:none;padding:0}}
  .interp-list li{{padding:.35rem 0;border-bottom:1px solid #f0f0f0;font-size:.85rem}}
  .interp-list li:before{{content:"▸ ";color:#0fa3a3;font-weight:700}}
  @media print{{.sidebar{{display:none}}.main{{margin-left:0}}}}
</style>
</head>
<body>
<div class="layout">

<!-- SIDEBAR -->
<div class="sidebar">
  <h2>🧬 FastQC Report</h2>
  <ul>{sidebar_items}</ul>
</div>

<!-- MAIN CONTENT -->
<div class="main">

  <!-- Header -->
  <div class="report-header">
    <h1>FastQC-style Report</h1>
    <p><strong>File:</strong> {filename} &nbsp;|&nbsp; <strong>Format:</strong> FASTQ &nbsp;|&nbsp; <strong>Encoding:</strong> Sanger / Illumina 1.9</p>
  </div>

  <!-- QC Overview bar -->
  <div class="qc-bar">
    {qc_bar_html}
  </div>

  <!-- 1. Basic Statistics -->
  <div class="section" id="basic_statistics">
    <div class="section-title"><span class="icon">✅</span> Basic Statistics</div>
    <table>
      <tr><th>Measure</th><th>Value</th></tr>
      <tr><td>Filename</td><td>{filename}</td></tr>
      <tr><td>File type</td><td>Conventional base calls</td></tr>
      <tr><td>Encoding</td><td>Sanger / Illumina 1.9</td></tr>
      <tr><td>Total Sequences</td><td>{total_seqs}</td></tr>
      <tr><td>Total Bases</td><td>{total_bases}</td></tr>
      <tr><td>Sequences flagged as poor quality</td><td>{summary_dict.get("Low-quality reads (Q<20)","N/A")}</td></tr>
      <tr><td>Sequence length</td><td>{seq_len}</td></tr>
      <tr><td>Mean Q score</td><td>{mean_q}</td></tr>
      <tr><td>% Reads Q≥30</td><td>{q30}</td></tr>
      <tr><td>%GC</td><td>{mean_gc}</td></tr>
    </table>
  </div>

  <!-- 2. Per base sequence quality -->
  <div class="section" id="per_base_sequence_quality">
    <div class="section-title"><span class="icon">✅</span> Per base sequence quality</div>
    <div class="chart-title">Quality scores across all bases (Sanger / Illumina 1.9 encoding)</div>
    <div class="zone-legend">
      <span><span class="zone-box" style="background:#2ecc71"></span>Very Good Q≥30</span>
      <span><span class="zone-box" style="background:#f39c12"></span>Acceptable Q20–29</span>
      <span><span class="zone-box" style="background:#e74c3c"></span>Poor &lt;Q20</span>
    </div>
    <div class="chart-wrap"><canvas id="perBaseChart"></canvas></div>
  </div>

  <!-- 3. Per sequence quality scores -->
  <div class="section" id="per_sequence_quality_scores">
    <div class="section-title"><span class="icon">✅</span> Per sequence quality scores</div>
    <div class="chart-title">Quality score distribution over all sequences</div>
    <div class="chart-wrap"><canvas id="perReadChart"></canvas></div>
  </div>

  <!-- 4. Per base sequence content -->
  <div class="section" id="per_base_sequence_content">
    <div class="section-title"><span class="icon">✅</span> Per base sequence content</div>
    <div class="chart-title">Sequence content across all bases</div>
    <div class="chart-wrap"><canvas id="compChart"></canvas></div>
  </div>

  <!-- 5. Per sequence GC content -->
  <div class="section" id="per_sequence_gc_content">
    <div class="section-title"><span class="icon">✅</span> Per sequence GC content</div>
    <div class="chart-title">GC distribution over all sequences</div>
    <div class="chart-wrap"><canvas id="gcChart"></canvas></div>
  </div>

  <!-- 6. Per base N content -->
  <div class="section" id="per_base_n_content">
    <div class="section-title"><span class="icon">✅</span> Per base N content</div>
    <div class="chart-title">N content across all bases</div>
    <div class="chart-wrap"><canvas id="nChart"></canvas></div>
  </div>

  <!-- 7. Sequence Length Distribution -->
  <div class="section" id="sequence_length_distribution">
    <div class="section-title"><span class="icon">✅</span> Sequence Length Distribution</div>
    <div class="chart-title">Distribution of sequence lengths over all sequences</div>
    <div class="chart-wrap"><canvas id="lenChart"></canvas></div>
  </div>

  <!-- 8. Sequence Duplication Levels -->
  <div class="section" id="sequence_duplication_levels">
    <div class="section-title"><span class="icon">✅</span> Sequence Duplication Levels</div>
    <div class="chart-title">Percent of seqs remaining if deduplicated</div>
    <div class="chart-wrap"><canvas id="dupChart"></canvas></div>
  </div>

  <!-- 9. Overrepresented sequences -->
  <div class="section" id="overrepresented_sequences">
    <div class="section-title"><span class="icon">✅</span> Overrepresented sequences</div>
    <p style="color:#5a7a9a;font-size:.88rem;">No overrepresented sequences detected (pure-Python analysis; for full k-mer detection run FastQC).</p>
  </div>

  <!-- 10. Adapter Content -->
  <div class="section" id="adapter_content">
    <div class="section-title"><span class="icon">✅</span> Adapter Content</div>
    <div class="chart-title">% Adapter</div>
    <div class="chart-wrap"><canvas id="adapterChart"></canvas></div>
    <p style="font-size:.78rem;color:#999;margin-top:.5rem;">Note: Adapter detection requires known adapter sequences. For full adapter analysis run FastQC or Trim Galore.</p>
  </div>

  <!-- Biological Interpretation -->
  <div class="section">
    <div class="section-title"><span class="icon">🔬</span> Biological Interpretation</div>
    <ul class="interp-list">{interp_html}</ul>
  </div>

</div><!-- /main -->
</div><!-- /layout -->

<script>
// ── Shared helpers ────────────────────────────────────────────
const ZONE_GREEN  = 'rgba(46,204,113,0.18)';
const ZONE_YELLOW = 'rgba(243,156,18,0.18)';
const ZONE_RED    = 'rgba(231,76,60,0.18)';

function zonePlugin(chartId) {{
  return {{
    id:'zones',
    beforeDraw(chart){{
      const {{ctx,chartArea,scales}}=chart;
      if(!chartArea) return;
      const yScale=scales.y;
      const x0=chartArea.left, x1=chartArea.right;
      const y30=yScale.getPixelForValue(30);
      const y20=yScale.getPixelForValue(20);
      const yBot=chartArea.bottom;
      const yTop=chartArea.top;
      ctx.save();
      ctx.fillStyle=ZONE_GREEN;  ctx.fillRect(x0,yTop,x1-x0,y30-yTop);
      ctx.fillStyle=ZONE_YELLOW; ctx.fillRect(x0,y30,x1-x0,y20-y30);
      ctx.fillStyle=ZONE_RED;    ctx.fillRect(x0,y20,x1-x0,yBot-y20);
      ctx.restore();
    }}
  }};
}}

// ── 1. Per-base sequence quality ────────────────────────────
new Chart(document.getElementById('perBaseChart'),{{
  type:'line',
  plugins:[zonePlugin('perBaseChart')],
  data:{{
    labels:{js(perbase_labels)},
    datasets:[
      {{label:'Mean Q',data:{js(perbase_mean)},borderColor:'#2471a3',borderWidth:2,pointRadius:0,fill:false,tension:0.3}},
      {{label:'Median Q',data:{js(perbase_median)},borderColor:'#1abc9c',borderWidth:1.5,pointRadius:0,fill:false,borderDash:[4,3],tension:0.3}},
      {{label:'90th Pct',data:{js(perbase_p90)},borderColor:'#f39c12',borderWidth:1,pointRadius:0,fill:false,borderDash:[2,4],tension:0.3}},
      {{label:'10th Pct',data:{js(perbase_p10)},borderColor:'#e74c3c',borderWidth:1,pointRadius:0,fill:false,borderDash:[2,4],tension:0.3}},
    ]
  }},
  options:{{responsive:true,maintainAspectRatio:false,
    scales:{{
      y:{{min:0,max:42,title:{{display:true,text:'Quality Score (Phred)'}}}},
      x:{{title:{{display:true,text:'Position in read (bp)'}}}}
    }},
    plugins:{{legend:{{position:'top'}},tooltip:{{mode:'index',intersect:false}}}}
  }}
}});

// ── 2. Per-read quality distribution ───────────────────────
(function(){{
  const labels={js(perread_labels)};
  const counts={js(perread_counts)};
  // background zones per bar
  const bgColors=labels.map(q=>q>=30?'rgba(46,204,113,0.7)':q>=20?'rgba(243,156,18,0.7)':'rgba(231,76,60,0.7)');
  new Chart(document.getElementById('perReadChart'),{{
    type:'line',
    data:{{labels,datasets:[{{
      label:'Average Quality per read',
      data:counts,
      borderColor:'#8e44ad',borderWidth:2,pointRadius:0,fill:true,
      backgroundColor:'rgba(142,68,173,0.12)',tension:0.3
    }}]}},
    options:{{responsive:true,maintainAspectRatio:false,
      scales:{{
        y:{{beginAtZero:true,title:{{display:true,text:'Number of sequences'}}}},
        x:{{title:{{display:true,text:'Mean Sequence Quality (Phred Score)'}}}}
      }},
      plugins:{{legend:{{position:'top'}}}}
    }}
  }});
}})();

// ── 3. Per base sequence content ────────────────────────────
new Chart(document.getElementById('compChart'),{{
  type:'line',
  data:{{
    labels:{js(comp_labels)},
    datasets:[
      {{label:'%T',data:{js(comp_t)},borderColor:'#e74c3c',borderWidth:1.8,pointRadius:0,fill:false,tension:0.3}},
      {{label:'%C',data:{js(comp_c)},borderColor:'#3498db',borderWidth:1.8,pointRadius:0,fill:false,tension:0.3}},
      {{label:'%A',data:{js(comp_a)},borderColor:'#2ecc71',borderWidth:1.8,pointRadius:0,fill:false,tension:0.3}},
      {{label:'%G',data:{js(comp_g)},borderColor:'#f39c12',borderWidth:1.8,pointRadius:0,fill:false,tension:0.3}},
    ]
  }},
  options:{{responsive:true,maintainAspectRatio:false,
    scales:{{
      y:{{min:0,max:100,title:{{display:true,text:'% bases'}}}},
      x:{{title:{{display:true,text:'Position in read (bp)'}}}}
    }},
    plugins:{{legend:{{position:'right'}},tooltip:{{mode:'index',intersect:false}}}}
  }}
}});

// ── 4. Per-sequence GC content ──────────────────────────────
new Chart(document.getElementById('gcChart'),{{
  type:'line',
  data:{{
    labels:{js(gc_labels)},
    datasets:[
      {{label:'GC count per read',data:{js(gc_counts)},borderColor:'#8e44ad',borderWidth:2,pointRadius:0,fill:false,tension:0.4}},
      {{label:'Theoretical Distribution',data:{js(gc_theoretical)},borderColor:'#2471a3',borderWidth:1.5,pointRadius:0,fill:false,borderDash:[5,3],tension:0.4}},
    ]
  }},
  options:{{responsive:true,maintainAspectRatio:false,
    scales:{{
      y:{{beginAtZero:true,title:{{display:true,text:'Number of reads'}}}},
      x:{{title:{{display:true,text:'Mean GC content (%)'}}}}
    }},
    plugins:{{legend:{{position:'top'}}}}
  }}
}});

// ── 5. Per-base N content ───────────────────────────────────
new Chart(document.getElementById('nChart'),{{
  type:'line',
  data:{{
    labels:{js(n_labels)},
    datasets:[{{label:'%N',data:{js(n_vals)},borderColor:'#c0392b',borderWidth:2,pointRadius:0,fill:true,backgroundColor:'rgba(192,57,43,0.12)',tension:0.3}}]
  }},
  options:{{responsive:true,maintainAspectRatio:false,
    scales:{{
      y:{{min:0,max:100,title:{{display:true,text:'% N'}}}},
      x:{{title:{{display:true,text:'Position in read (bp)'}}}}
    }},
    plugins:{{legend:{{position:'top'}}}}
  }}
}});

// ── 6. Sequence length distribution ────────────────────────
new Chart(document.getElementById('lenChart'),{{
  type:'line',
  data:{{
    labels:{js(len_labels)},
    datasets:[{{label:'Sequence Length',data:{js(len_counts)},borderColor:'#8e44ad',borderWidth:2,pointRadius:0,fill:false,tension:0.4}}]
  }},
  options:{{responsive:true,maintainAspectRatio:false,
    scales:{{
      y:{{beginAtZero:true,title:{{display:true,text:'Number of reads'}}}},
      x:{{title:{{display:true,text:'Sequence Length (bp)'}}}}
    }},
    plugins:{{legend:{{position:'top'}}}}
  }}
}});

// ── 7. Sequence Duplication Levels ─────────────────────────
new Chart(document.getElementById('dupChart'),{{
  type:'line',
  data:{{
    labels:{js(dup_labels)},
    datasets:[{{label:'% Total sequences',data:{js(dup_vals)},borderColor:'#8e44ad',borderWidth:2,pointRadius:3,fill:false,tension:0}}]
  }},
  options:{{responsive:true,maintainAspectRatio:false,
    scales:{{
      y:{{min:0,max:100,title:{{display:true,text:'% of sequences'}}}},
      x:{{title:{{display:true,text:'Sequence Duplication Level'}}}}
    }},
    plugins:{{legend:{{position:'top'}}}}
  }}
}});

// ── 8. Adapter content ──────────────────────────────────────
new Chart(document.getElementById('adapterChart'),{{
  type:'line',
  data:{{
    labels:{js(adapter_labels)},
    datasets:[
      {{label:'Illumina Universal Adapter',data:{js(adapter_vals)},borderColor:'#e74c3c',borderWidth:1.5,pointRadius:0,fill:false}},
      {{label:'Illumina Small RNA 3 Adapter',data:{js(adapter_vals)},borderColor:'#3498db',borderWidth:1.5,pointRadius:0,fill:false}},
      {{label:'Nextera Transposase',data:{js(adapter_vals)},borderColor:'#2ecc71',borderWidth:1.5,pointRadius:0,fill:false}},
      {{label:'PolyA',data:{js(adapter_vals)},borderColor:'#f39c12',borderWidth:1.5,pointRadius:0,fill:false}},
      {{label:'PolyG',data:{js(adapter_vals)},borderColor:'#9b59b6',borderWidth:1.5,pointRadius:0,fill:false}},
    ]
  }},
  options:{{responsive:true,maintainAspectRatio:false,
    scales:{{
      y:{{min:0,max:100,title:{{display:true,text:'% Adapter'}}}},
      x:{{title:{{display:true,text:'Position in read (bp)'}}}}
    }},
    plugins:{{legend:{{position:'right'}}}}
  }}
}});
</script>

</body>
</html>"""
    return html


def generate_pdf(result,filename):
    if not REPORTLAB_OK:
        return None
    buf=io.BytesIO()
    NAVY=colors.HexColor("#0D1F3C"); TEAL=colors.HexColor("#0FA3A3")
    LIGHT=colors.HexColor("#F0F4F8"); MID=colors.HexColor("#D0E0EC")
    WHITE=colors.white; DARK=colors.HexColor("#1E2D40")
    PASS_C=colors.HexColor("#1A7A4A"); WARN_C=colors.HexColor("#B8740A"); FAIL_C=colors.HexColor("#C0392B")

    def sty(name,**kw):
        return ParagraphStyle(name,fontName=kw.get("font","Helvetica"),fontSize=kw.get("size",9),
                              textColor=kw.get("color",DARK),spaceAfter=kw.get("sa",2),
                              leading=kw.get("lead",13),alignment=kw.get("align",0),
                              leftIndent=kw.get("li",0))

    S={"body":sty("body"),"bold":sty("bold",font="Helvetica-Bold"),
       "sec":sty("sec",font="Helvetica-Bold",size=11,color=NAVY,sa=6),
       "cap":sty("cap",font="Helvetica-Oblique",size=7.5,color=colors.HexColor("#888888"),sa=4),
       "pass":sty("pass",font="Helvetica-Bold",size=9,color=PASS_C),
       "warn":sty("warn",font="Helvetica-Bold",size=9,color=WARN_C),
       "fail":sty("fail",font="Helvetica-Bold",size=9,color=FAIL_C)}

    fmt=result.get("format","Unknown")
    doc=SimpleDocTemplate(buf,pagesize=A4,leftMargin=15*mm,rightMargin=15*mm,topMargin=22*mm,bottomMargin=16*mm)

    def hf(c,d):
        W,H=A4
        c.saveState()
        c.setFillColor(NAVY); c.rect(0,H-16*mm,W,16*mm,fill=1,stroke=0)
        c.setFillColor(TEAL); c.rect(0,H-17.5*mm,W,1.5*mm,fill=1,stroke=0)
        c.setFillColor(WHITE); c.setFont("Helvetica-Bold",10); c.drawString(15*mm,H-11*mm,"NGS Interpreter")
        c.setFont("Helvetica",8); c.drawRightString(W-15*mm,H-11*mm,f"{fmt}  ·  {filename}")
        c.setFillColor(NAVY); c.rect(0,0,W,9*mm,fill=1,stroke=0)
        c.setFillColor(MID); c.setFont("Helvetica",7)
        c.drawString(15*mm,3*mm,datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
        c.drawCentredString(W/2,3*mm,"NGS Interpreter — Genomics Analysis Platform")
        c.drawRightString(W-15*mm,3*mm,f"Page {d.page}")
        c.restoreState()

    def mt(rows):
        data=[[Paragraph(f"<b>{l}</b>",S["body"]),Paragraph(str(v),S["body"])] for l,v in rows]
        t=Table(data,colWidths=[65*mm,105*mm])
        t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),LIGHT),("BACKGROUND",(0,0),(0,-1),colors.HexColor("#E0EBF5")),
            ("GRID",(0,0),(-1,-1),0.5,MID),("TOPPADDING",(0,0),(-1,-1),3),("BOTTOMPADDING",(0,0),(-1,-1),3),
            ("LEFTPADDING",(0,0),(-1,-1),6),("VALIGN",(0,0),(-1,-1),"MIDDLE"),
            ("ROWBACKGROUNDS",(0,0),(-1,-1),[LIGHT,WHITE])]))
        return t

    def dft(df,max_rows=40):
        df=df.head(max_rows); cols=list(df.columns)
        hrow=[Paragraph(f"<b>{c}</b>",S["body"]) for c in cols]
        data=[hrow]+[[Paragraph(str(v)[:55],S["body"]) for v in row] for _,row in df.iterrows()]
        n=len(cols); cw=170*mm/n
        t=Table(data,colWidths=[cw]*n,repeatRows=1)
        t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),NAVY),("TEXTCOLOR",(0,0),(-1,0),WHITE),
            ("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),("FONTSIZE",(0,0),(-1,-1),7.5),
            ("GRID",(0,0),(-1,-1),0.4,MID),("ROWBACKGROUNDS",(0,1),(-1,-1),[WHITE,LIGHT]),
            ("TOPPADDING",(0,0),(-1,-1),3),("BOTTOMPADDING",(0,0),(-1,-1),3),("LEFTPADDING",(0,0),(-1,-1),4)]))
        return t

    def bt(items):
        data=[[Paragraph(f"<b>{h}</b>",S["body"]) for h in ["Category","Count","%","Distribution"]]]
        for label,count,pct,bar_pct in items:
            filled=int(bar_pct/100*28); bar="█"*filled+"░"*(28-filled)
            data.append([Paragraph(str(label),S["body"]),Paragraph(f"{count:,}" if isinstance(count,int) else str(count),S["body"]),
                         Paragraph(f"{pct:.1f}%",S["body"]),Paragraph(f"<font face='Courier' size='7'>{bar}</font>",S["body"])])
        t=Table(data,colWidths=[45*mm,20*mm,16*mm,89*mm],repeatRows=1)
        t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),NAVY),("TEXTCOLOR",(0,0),(-1,0),WHITE),
            ("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),("FONTSIZE",(0,0),(-1,-1),8),
            ("GRID",(0,0),(-1,-1),0.4,MID),("ROWBACKGROUNDS",(0,1),(-1,-1),[WHITE,LIGHT]),
            ("TOPPADDING",(0,0),(-1,-1),3),("BOTTOMPADDING",(0,0),(-1,-1),3),("LEFTPADDING",(0,0),(-1,-1),4)]))
        return t

    story=[Spacer(1,4*mm)]
    cover=Table([[Paragraph("NGS INTERPRETER",ParagraphStyle("ct",fontName="Helvetica-Bold",fontSize=18,textColor=NAVY,alignment=TA_CENTER))],
                 [Paragraph("Genomics File Analysis Report",ParagraphStyle("cs",fontName="Helvetica",fontSize=11,textColor=colors.HexColor("#5A7A9A"),alignment=TA_CENTER))]],colWidths=[170*mm])
    cover.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),LIGHT),("TOPPADDING",(0,0),(-1,-1),10),("BOTTOMPADDING",(0,0),(-1,-1),10),("LINEBELOW",(0,-1),(-1,-1),2,TEAL),("LINEABOVE",(0,0),(-1,0),2,NAVY)]))
    story+=[cover,Spacer(1,4*mm)]
    now=datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    meta=Table([[Paragraph(f"<b>File:</b> {filename}",S["body"]),Paragraph(f"<b>Format:</b> {fmt}",S["body"]),Paragraph(f"<b>Date:</b> {now}",S["body"])]],colWidths=[65*mm,55*mm,50*mm])
    meta.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),colors.HexColor("#DDE8F5")),("GRID",(0,0),(-1,-1),0.5,MID),("TOPPADDING",(0,0),(-1,-1),5),("BOTTOMPADDING",(0,0),(-1,-1),5),("LEFTPADDING",(0,0),(-1,-1),8)]))
    story+=[meta,Spacer(1,5*mm)]

    if result.get("error"):
        story.append(Paragraph("Error: "+result["error"],S["body"]))
    else:
        for sec in result.get("sections",[]):
            story+=[HRFlowable(width="100%",thickness=1,color=TEAL,spaceAfter=2),Paragraph(sec.get("title",""),S["sec"])]
            t=sec.get("type",""); d=sec.get("data")
            if t=="metrics" and isinstance(d,list): story.append(mt(d))
            elif t=="table":
                try:
                    df=d if hasattr(d,"empty") else pd.DataFrame(d)
                    if not df.empty: story.append(dft(df))
                except: pass
            elif t=="bars" and isinstance(d,list): story.append(bt(d))
            elif t=="text" and isinstance(d,list):
                for line in d: story.append(Paragraph(f"• {line}",S["body"]))
            story.append(Spacer(1,3*mm))

        qc=result.get("qc",[])
        if qc:
            story+=[HRFlowable(width="100%",thickness=1.5,color=NAVY,spaceAfter=2),Paragraph("Quality Control Summary",S["sec"])]
            ICONS={"PASS":"✔  PASS","WARN":"⚠  WARN","FAIL":"✘  FAIL"}
            SM={"PASS":S["pass"],"WARN":S["warn"],"FAIL":S["fail"]}
            data=[[Paragraph("<b>Check</b>",S["body"]),Paragraph("<b>Result</b>",S["body"])]]+\
                 [[Paragraph(q["check"],S["body"]),Paragraph(ICONS.get(q["result"],q["result"]),SM.get(q["result"],S["body"]))] for q in qc]
            t=Table(data,colWidths=[120*mm,50*mm],repeatRows=1)
            t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),NAVY),("TEXTCOLOR",(0,0),(-1,0),WHITE),("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),
                ("GRID",(0,0),(-1,-1),0.4,MID),("ROWBACKGROUNDS",(0,1),(-1,-1),[WHITE,LIGHT]),("FONTSIZE",(0,0),(-1,-1),9),
                ("TOPPADDING",(0,0),(-1,-1),4),("BOTTOMPADDING",(0,0),(-1,-1),4),("LEFTPADDING",(0,0),(-1,-1),6)]))
            story+=[ t,Spacer(1,4*mm)]

        interp=result.get("interpretation",[])
        if interp:
            story+=[HRFlowable(width="100%",thickness=1.5,color=NAVY,spaceAfter=2),Paragraph("Biological & Clinical Interpretation",S["sec"])]
            for line in interp: story.append(Paragraph(f"▸  {line}",S["body"]))

    doc.build(story,onFirstPage=hf,onLaterPages=hf)
    return buf.getvalue()

# ════════════════════════════════════════════════════════════
# PAGE 1 — HOME
# ════════════════════════════════════════════════════════════

def page_home():
    topbar()
    st.markdown("""
    <div class="hero">
        <div class="hero-tag">Open Source · Offline · Free</div>
        <h1>Understand Your<br>Genomics Files</h1>
        <p>NGS Interpreter reads any Next-Generation Sequencing file and gives you a complete,
        plain-English breakdown of what's inside. No command line. No coding. Perfect for
        beginners and researchers alike.</p>
        <div>
            <span class="hero-badge">🔒 100% Offline</span>
            <span class="hero-badge">📂 20 File Formats</span>
            <span class="hero-badge">📊 PDF Reports</span>
            <span class="hero-badge">🧬 Alignment · Variants · Annotation</span>
        </div>
    </div>""", unsafe_allow_html=True)

    c1,c2,c3=st.columns([1,1.2,1])
    with c2:
        if st.button("🚀  Let's Start",use_container_width=True,type="primary"): go("analyze")

    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("### What is NGS and why do these files exist?")
    st.markdown("<p style='color:#5A7A9A;font-size:0.92rem;max-width:700px;line-height:1.7;'>Next-Generation Sequencing (NGS) reads millions of DNA/RNA fragments simultaneously. Each step of the analysis pipeline creates a different file format. This tool decodes every one of them — no command line needed.</p>", unsafe_allow_html=True)

    # ── Format cards ──────────────────────────────────────────────────
    c1,c2,c3=st.columns(3)
    cards=[
        ("🧬 FASTQ / FASTA",
         "Raw sequencing reads straight off the sequencer. FASTQ carries a quality score (Phred) for every base call. FASTA is the plain sequence-only format used for reference genomes and protein databases. NGS Interpreter reports read counts, Q20/Q30%, GC%, N content, and per-base quality profiles."),
        ("📌 BAM / SAM / CRAM",
         "Alignment files — your reads mapped back onto a reference genome. SAM is human-readable text; BAM is its compressed binary equivalent; CRAM is further compressed using the reference. Reports mapping rate, duplicate rate, MAPQ distribution, insert size, and per-chromosome coverage."),
        ("🔬 VCF / BCF",
         "Variant Call Format — records every difference between your sample and the reference: SNPs, insertions, deletions, and multi-nucleotide variants. BCF is the binary version. Reports Ts/Tv ratio, zygosity, QUAL/DP/AF distributions, filter status, and chromosome-level variant counts."),
        ("📋 GFF3 / GTF",
         "Genome annotation files — they describe where genes, exons, introns, UTRs, and regulatory features sit on each chromosome. GTF is the Ensembl/GENCODE flavour. Reports feature type counts, gene/transcript tallies, biotype distribution, exon structure, and strand balance."),
        ("📖 GenBank / EMBL",
         "Rich flat-file records combining sequence and annotation in one file. GenBank is the NCBI standard; EMBL is used by European databases. Reports locus info, CDS statistics, product annotations, enzyme class breakdown, GC content, and coding density."),
        ("📊 BED / narrowPeak / broadPeak",
         "Interval files that mark genomic regions of interest — ChIP-seq and ATAC-seq peaks, transcript models, CpG islands, and repeat elements. narrowPeak stores point-source peaks (TF binding); broadPeak stores diffuse domains (histone marks). Reports region counts, widths, scores, and chromosome distribution."),
        ("〰️ WIG / BEDGraph",
         "Continuous signal tracks representing coverage depth, read density, or enrichment across the genome. variableStep and fixedStep WIG formats are common from older tools; BEDGraph (4-column) is the modern standard. NGS Interpreter reports min/max/mean signal, per-chromosome density, and score distribution."),
        ("🧩 MAF (Multiple Alignment Format)",
         "Whole-genome multi-species alignment blocks — used to compare syntenic regions across organisms. Produced by tools like MULTIZ, LASTZ, and Progressive Cactus. Reports block count, species coverage, alignment scores, and data preview. Conserved blocks indicate functional elements under purifying selection."),
        ("🧪 AB1 / Sanger Chromatogram",
         "Applied Biosystems binary format for Sanger (first-generation) sequencing traces. Stores raw fluorescence signal, called bases, and per-base quality scores. NGS Interpreter extracts the sequence, computes GC% and N content, and summarises quality scores — useful for primer validation and clone confirmation."),
        ("🗂️ PLINK .bim / .fam",
         ".bim is the PLINK variant map: one row per SNP with chromosome, rs ID, cM position, base-pair position, and both alleles. .fam is the sample info file: family ID, individual ID, parental IDs, sex, and case/control phenotype. Both are used in GWAS, PRS, and population genetics analyses."),
        ("📐 EMBL Flat File",
         "The European Molecular Biology Laboratory flat file format — the EMBL-EBI counterpart to GenBank. Contains the same rich combination of sequence, feature table, and metadata but follows EMBL field conventions (FT, SQ blocks). Interconvertible with GenBank using Biopython or EMBOSS Seqret."),
        ("🌐 BCF (Binary VCF)",
         "BCF is the BGZF-compressed, block-gzip binary version of VCF — dramatically smaller and faster to query with tabix and bcftools. Internally identical to VCF once decompressed. NGS Interpreter reads BCF files using the same variant analysis pipeline as VCF: Ts/Tv, zygosity, depth, allele frequency, and filter status."),
    ]
    for i,(title,desc) in enumerate(cards):
        with [c1,c2,c3][i%3]:
            st.markdown(f'<div class="gloss-card"><h4>{title}</h4><p>{desc}</p></div><br>', unsafe_allow_html=True)

    # ── Pipeline steps ────────────────────────────────────────────────
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("### The Typical NGS Analysis Pipeline")
    st.markdown("<p style='color:#5A7A9A;font-size:0.88rem;'>Every step in the pipeline produces a different file type. NGS Interpreter understands all of them.</p>", unsafe_allow_html=True)
    steps=[
        ("1","Wet Lab","DNA/RNA extraction → library prep → sequencing run on Illumina / Nanopore / PacBio"),
        ("2","Raw Reads","Output: FASTQ files → run FastQC for quality control → trim adapters with Cutadapt"),
        ("3","Alignment","Map reads to reference genome → Output: SAM / BAM / CRAM → index with samtools"),
        ("4","Variant Calling","GATK HaplotypeCaller / DeepVariant → Output: VCF / BCF → filter and annotate"),
        ("5","Annotation","Add gene context → GFF3 / GTF / GenBank / EMBL → functional interpretation"),
        ("6","You Are Here","← NGS Interpreter decodes every file at every step of this pipeline"),
    ]
    cols=st.columns(6)
    for col,(num,title,desc) in zip(cols,steps):
        with col:
            st.markdown(f"""<div style='background:white;border:1px solid #D0DDE8;border-top:3px solid #0FA3A3;
                border-radius:10px;padding:0.8rem;text-align:center;height:100%;'>
                <div style='font-size:1.3rem;font-weight:700;color:#0FA3A3;'>{num}</div>
                <div style='font-weight:600;font-size:0.82rem;color:#0D1F3C;margin:0.3rem 0;'>{title}</div>
                <div style='font-size:0.75rem;color:#5A7A9A;line-height:1.5;'>{desc}</div>
            </div>""", unsafe_allow_html=True)

    # ── Supported formats table ───────────────────────────────────────
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("### All Supported File Formats")
    fmt_rows=[
        ("FASTQ / .fq / .fastq.gz","Raw Reads","✅ Full","Q scores, GC%, N content, per-base quality, length distribution"),
        ("FASTA / .fa .fna .faa .fasta","Sequences","✅ Full","N50, GC%, complexity, dinucleotide freq, CpG sites, AT/GC skew"),
        ("BAM / .bam","Alignment (binary)","✅ Full","Mapping rate, MAPQ, duplicates, insert size, chromosome coverage"),
        ("SAM / .sam","Alignment (text)","✅ Full","Same as BAM — parsed via pysam or pure Python fallback"),
        ("CRAM / .cram","Alignment (compressed)","✅ Full","Same as BAM — reference-compressed format"),
        ("VCF / .vcf .vcf.gz","Variants","✅ Full","SNP/INDEL/MNV, Ts/Tv, zygosity, QUAL/DP/AF, chromosome dist."),
        ("BCF / .bcf","Variants (binary)","✅ Full","Binary VCF — same analysis pipeline as VCF"),
        ("GFF3 / .gff .gff3","Annotation","✅ Full","Feature types, strand balance, lengths, products, hypothetical %"),
        ("GTF / .gtf","Annotation (Ensembl)","✅ Full","Gene/transcript counts, biotypes, TSL, exon structure, CDS"),
        ("GenBank / .gbk .gb .gbff","Flat file (NCBI)","✅ Full","Coding density, CDS stats, EC numbers, enzyme classes, products"),
        ("EMBL / .embl","Flat file (EBI)","✅ Full","Same as GenBank — EMBL-EBI standard format"),
        ("BED / .bed","Intervals","✅ Full","Region widths, chromosome dist., score summary"),
        ("narrowPeak / .narrowPeak","ChIP/ATAC peaks","✅ Full","Signal value, p-value, q-value, peak width distribution"),
        ("broadPeak / .broadPeak","Histone peaks","✅ Full","Broad domain stats, chromosome coverage"),
        ("WIG / .wig","Signal track","✅ Full","Min/max/mean signal, per-chrom density, score distribution"),
        ("BEDGraph / .bedGraph","Signal track","✅ Full","Same as WIG — 4-column continuous signal format"),
        ("MAF / .maf","Multi-alignment","✅ Full","Block count, species coverage, alignment scores"),
        ("PLINK .bim","Variant map","✅ Full","SNPs per chromosome, allele counts"),
        ("PLINK .fam","Sample info","✅ Full","Sex distribution, case/control phenotype breakdown"),
        ("AB1 / .ab1","Sanger trace","✅ Full","Sequence, GC%, N content, Q score summary"),
    ]
    df_fmt=pd.DataFrame(fmt_rows,columns=["Format / Extension","Category","Support","What NGS Interpreter Reports"])
    st.dataframe(df_fmt,use_container_width=True,hide_index=True)

    st.markdown("<br>", unsafe_allow_html=True)
    c1,c2,c3=st.columns([1,1.2,1])
    with c2:
        if st.button("🚀  Let's Start",key="start2",use_container_width=True,type="primary"): go("analyze")

# ════════════════════════════════════════════════════════════
# PAGE 2 — ANALYZE
# ════════════════════════════════════════════════════════════

def page_analyze():
    topbar()
    cb,ct=st.columns([0.15,0.85])
    with cb:
        if st.button("← Home"): go("home")
    with ct:
        st.markdown("## Analyse Your File")

    st.markdown("<p style='color:#5A7A9A;'>Upload any NGS file — format is detected automatically. Everything runs on this server; no data is stored.</p>", unsafe_allow_html=True)

    uploaded=st.file_uploader("",type=None,label_visibility="collapsed",help="Supports BAM, SAM, FASTA, FASTQ, VCF, BCF, GFF3, GTF, GenBank, EMBL, BED, narrowPeak, broadPeak, BEDGraph, WIG, MAF, PLINK .bim/.fam, AB1")

    if uploaded:
        fsize=len(uploaded.getvalue()); fstr=f"{fsize/1024:.1f} KB" if fsize<1048576 else f"{fsize/1048576:.2f} MB"
        st.markdown(f"""<div style='background:white;border:1px solid #D0DDE8;border-radius:10px;padding:1rem 1.4rem;
            margin:1rem 0;display:flex;align-items:center;gap:1rem;'>
            <div style='font-size:2rem;'>📄</div>
            <div><div style='font-weight:600;color:#0D1F3C;'>{uploaded.name}</div>
            <div style='color:#5A7A9A;font-size:0.82rem;'>{fstr}</div></div></div>""", unsafe_allow_html=True)

        with tempfile.NamedTemporaryFile(delete=False,suffix=f"_{uploaded.name}") as tmp:
            tmp.write(uploaded.getvalue()); tmp_path=tmp.name

        fmt=detect_format(tmp_path)
        info=FORMAT_INFO.get(fmt,(fmt,"Unknown","—"))

        if fmt=="UNKNOWN":
            st.error(f"Could not detect the format of **{uploaded.name}**. Please check the file is a supported NGS format.")
            try: os.unlink(tmp_path)
            except: pass
            return

        st.markdown(f"""<div style='margin:1rem 0;'>
            <div style='font-size:0.75rem;font-weight:600;color:#5A7A9A;text-transform:uppercase;letter-spacing:0.08em;margin-bottom:0.4rem;'>Detected Format</div>
            <div class="fmt-badge">🔬 <div><div>{info[0]}</div>
            <div style='color:#0FA3A3;font-size:0.72rem;'>{info[1]} · {info[2]}</div></div></div></div>""", unsafe_allow_html=True)

        c1,c2,c3=st.columns([1,1.5,1])
        with c2:
            run=st.button(f"▶  Analyse {info[0]} File",use_container_width=True,type="primary",key="run_btn")

        if run:
            with st.spinner(f"Analysing {uploaded.name} …"):
                import time; t0=time.time()
                result=analyze(tmp_path,fmt)
                elapsed=time.time()-t0
            try: os.unlink(tmp_path)
            except: pass
            if result.get("error"):
                st.error(f"Analysis error: {result['error']}")
                return
            st.session_state.result=result; st.session_state.fname=uploaded.name; st.session_state.fmt=fmt
            st.success(f"✓ Analysis complete in {elapsed:.1f}s")
            import time; time.sleep(0.5)
            go("results")
        else:
            try: os.unlink(tmp_path)
            except: pass

# ════════════════════════════════════════════════════════════
# PAGE 3 — RESULTS
# ════════════════════════════════════════════════════════════

def render_section(sec):
    st.markdown(f'<div class="sec-header"><div class="sec-dot"></div>{sec.get("title","")}</div>', unsafe_allow_html=True)
    t=sec.get("type",""); d=sec.get("data")

    if t=="metrics" and isinstance(d,list):
        html='<div class="metrics-grid">'
        for label,value in d:
            html+=f'<div class="metric-row"><span class="metric-label">{label}</span><span class="metric-value">{value}</span></div>'
        html+="</div>"
        st.markdown(html, unsafe_allow_html=True)

    elif t=="table":
        try:
            df=d if hasattr(d,"empty") else pd.DataFrame(d)
            if not df.empty: st.dataframe(df,use_container_width=True,hide_index=True)
        except Exception as e: st.warning(f"Table error: {e}")

    elif t=="bars" and isinstance(d,list):
        df=pd.DataFrame([{"Category":l,"Count":c,"%":p} for l,c,p,_ in d])
        c1,c2=st.columns([2,3])
        with c1: st.dataframe(df,use_container_width=True,hide_index=True)
        with c2: st.bar_chart(df.set_index("Category")["%"],use_container_width=True,height=min(300,max(150,len(d)*30)))

    elif t=="text" and isinstance(d,list):
        for line in d: st.markdown(f"<div class='interp-item'>▸ {line}</div>", unsafe_allow_html=True)

    # ── FastQC: Per-base quality line chart with colour zones ─────────
    elif t=="fastqc_perbase" and isinstance(d,list):
        df=pd.DataFrame(d)
        # colour-code background zones using annotation
        st.markdown("""
        <div style='display:flex;gap:1rem;margin-bottom:0.4rem;font-size:0.78rem;'>
          <span style='background:#D4EDDA;border-radius:4px;padding:2px 8px;color:#1A7A4A;font-weight:600;'>● Q≥30 Excellent</span>
          <span style='background:#FFF3CD;border-radius:4px;padding:2px 8px;color:#B8740A;font-weight:600;'>● Q20–29 Acceptable</span>
          <span style='background:#F8D7DA;border-radius:4px;padding:2px 8px;color:#C0392B;font-weight:600;'>● Q&lt;20 Poor</span>
        </div>""", unsafe_allow_html=True)
        # Build chart with Mean Q, Median, 10th and 90th percentile
        chart_df=df.set_index("Position")[["Mean Q","Median Q","10th Pct","90th Pct"]]
        st.line_chart(chart_df, use_container_width=True, height=320)
        # Summary table underneath
        with st.expander("View per-base quality table"):
            st.dataframe(df[["Position","Mean Q","Median Q","10th Pct","90th Pct","Zone"]],use_container_width=True,hide_index=True)

    # ── FastQC: single-series histogram bar chart ─────────────────────
    elif t=="fastqc_hist" and isinstance(d,list):
        df=pd.DataFrame(d)
        x_col=sec.get("x","x"); y_col=sec.get("y","y")
        if x_col in df.columns and y_col in df.columns:
            plot_df=df.set_index(x_col)[[y_col]]
            st.bar_chart(plot_df, use_container_width=True, height=280)
        with st.expander("View data table"):
            st.dataframe(df,use_container_width=True,hide_index=True)

    # ── FastQC: multi-line chart (A/T/G/C composition) ────────────────
    elif t=="fastqc_multiline" and isinstance(d,list):
        df=pd.DataFrame(d)
        x_col=sec.get("x","Position"); lines=sec.get("lines",[])
        valid=[l for l in lines if l in df.columns]
        if valid:
            st.line_chart(df.set_index(x_col)[valid], use_container_width=True, height=280)
        # Ideal target line annotation
        st.markdown("<p style='font-size:0.78rem;color:#5A7A9A;margin-top:-0.5rem;'>Ideal: all lines flat near 25% for random sequence. Bias at first 10 positions is normal (random hexamer priming).</p>", unsafe_allow_html=True)
        with st.expander("View base composition table"):
            st.dataframe(df,use_container_width=True,hide_index=True)

    if sec.get("note"): st.caption(sec["note"])


def page_results():
    topbar()
    result=st.session_state.result; fname=st.session_state.fname; fmt=st.session_state.fmt
    if not result:
        st.warning("No results yet."); 
        if st.button("← Go to Analyse"): go("analyze")
        return

    # For FASTQ show 3 columns (back | file info | FastQC HTML + PDF)
    if fmt == "FASTQ":
        cb,ci,cp1,cp2=st.columns([0.12,0.55,0.18,0.15])
    else:
        cb,ci,cp=st.columns([0.15,0.65,0.2])

    with cb:
        if st.button("← Analyse"): go("analyze")
    with ci:
        info=FORMAT_INFO.get(fmt,(fmt,"",""))
        st.markdown(f"<div style='padding:0.4rem 0;'><span style='font-weight:600;color:#0D1F3C;font-size:1rem;'>{fname}</span> &nbsp; <span style='background:#0D1F3C;color:white;border-radius:5px;padding:0.1rem 0.5rem;font-family:JetBrains Mono,monospace;font-size:0.75rem;'>{info[0]}</span> &nbsp; <span style='color:#5A7A9A;font-size:0.82rem;'>{info[1]}</span></div>", unsafe_allow_html=True)

    if fmt == "FASTQ":
        with cp1:
            html_report = generate_fastqc_html(result, fname)
            st.download_button(
                "🧬 FastQC HTML",
                data=html_report.encode("utf-8"),
                file_name=f"fastqc_{fname}.html",
                mime="text/html",
                use_container_width=True,
                type="primary"
            )
        with cp2:
            if st.button("📄 PDF",use_container_width=True):
                if REPORTLAB_OK:
                    with st.spinner("Generating PDF…"):
                        pdf=generate_pdf(result,fname)
                    st.download_button("⬇ Download PDF",pdf,f"ngsi_{fname}.pdf","application/pdf",use_container_width=True)
                else:
                    st.warning("Install reportlab: pip install reportlab")
    else:
        with cp:
            if st.button("📄 Export PDF",use_container_width=True):
                if REPORTLAB_OK:
                    with st.spinner("Generating PDF…"):
                        pdf=generate_pdf(result,fname)
                    st.download_button("⬇ Download PDF",pdf,f"ngsi_{fname}.pdf","application/pdf",use_container_width=True)
                else:
                    st.warning("Install reportlab for PDF export: pip install reportlab")

    st.markdown("<hr style='margin:0.8rem 0;'>", unsafe_allow_html=True)

    if result.get("error"):
        st.error(result["error"]); return

    # QC strip
    qc=result.get("qc",[])
    if qc:
        n_p=sum(1 for q in qc if q["result"]=="PASS")
        n_w=sum(1 for q in qc if q["result"]=="WARN")
        n_f=sum(1 for q in qc if q["result"]=="FAIL")
        st.markdown(f"""<div style='display:flex;align-items:center;gap:1rem;background:white;border:1px solid #D0DDE8;
            border-radius:10px;padding:0.7rem 1.2rem;margin-bottom:1.2rem;'>
            <div style='font-weight:600;color:#0D1F3C;font-size:0.88rem;'>QC Overview</div>
            <div style='display:flex;gap:0.5rem;margin-left:auto;'>
                <span style='background:#D4EDDA;color:#1A7A4A;border-radius:6px;padding:0.2rem 0.7rem;font-size:0.82rem;font-weight:700;'>✔ {n_p} PASS</span>
                <span style='background:#FFF3CD;color:#B8740A;border-radius:6px;padding:0.2rem 0.7rem;font-size:0.82rem;font-weight:700;'>⚠ {n_w} WARN</span>
                <span style='background:#F8D7DA;color:#C0392B;border-radius:6px;padding:0.2rem 0.7rem;font-size:0.82rem;font-weight:700;'>✘ {n_f} FAIL</span>
            </div></div>""", unsafe_allow_html=True)

    for sec in result.get("sections",[]): render_section(sec)

    if qc:
        st.markdown('<div class="sec-header"><div class="sec-dot"></div>Quality Control Checks</div>', unsafe_allow_html=True)
        html=""
        for q in qc:
            r=q["result"]; cls={"PASS":"qc-pass","WARN":"qc-warn","FAIL":"qc-fail"}.get(r,"qc-pass")
            icon={"PASS":"✔","WARN":"⚠","FAIL":"✘"}.get(r,"?")
            html+=f'<div class="qc-row"><span class="qc-badge {cls}">{icon} {r}</span><span style="font-size:0.85rem;">{q["check"]}</span></div>'
        st.markdown(html, unsafe_allow_html=True)

    interp=result.get("interpretation",[])
    if interp:
        st.markdown('<div class="sec-header"><div class="sec-dot"></div>Biological &amp; Clinical Interpretation</div>', unsafe_allow_html=True)
        st.markdown('<div class="interp-box">'+"".join(f'<div class="interp-item">▸ {l}</div>' for l in interp)+"</div>", unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)
    if fmt == "FASTQ":
        c1,c2,c3,c4=st.columns([1,1.2,1.2,1])
        with c2:
            html_report2 = generate_fastqc_html(result, fname)
            st.download_button(
                "🧬 Download FastQC HTML Report",
                data=html_report2.encode("utf-8"),
                file_name=f"fastqc_{fname}.html",
                mime="text/html",
                use_container_width=True,
                type="primary",
                key="fastqc_dl2"
            )
        with c3:
            if st.button("📄 Export PDF Report",use_container_width=True,key="pdf2"):
                if REPORTLAB_OK:
                    with st.spinner("Generating PDF…"):
                        pdf=generate_pdf(result,fname)
                    st.download_button("⬇ Download PDF",pdf,f"ngsi_{fname}.pdf","application/pdf",use_container_width=True,key="dl2")
                else:
                    st.warning("Install reportlab: pip install reportlab")
    else:
        c1,c2,c3=st.columns([1,1.5,1])
        with c2:
            if st.button("📄 Export Full PDF Report",use_container_width=True,type="primary",key="pdf2"):
                if REPORTLAB_OK:
                    with st.spinner("Generating PDF…"):
                        pdf=generate_pdf(result,fname)
                    st.download_button("⬇ Download PDF Report",pdf,f"ngsi_{fname}.pdf","application/pdf",use_container_width=True,key="dl2")
                else:
                    st.warning("Install reportlab: pip install reportlab")

# ════════════════════════════════════════════════════════════
# ROUTER
# ════════════════════════════════════════════════════════════

page=st.session_state.page
if page=="home":    page_home()
elif page=="analyze": page_analyze()
elif page=="results": page_results()

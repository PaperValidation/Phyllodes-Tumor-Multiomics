import csv, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from _env import TABLES as TAB, DATASET_PATH

def cidx(header, name):
    for i,h in enumerate(header):
        if h.strip()==name: return i
    return -1

def sig_adjP(path):
    with open(path, newline='') as f:
        r=csv.reader(f); h=next(r); idx=cidx(h,"adj.P.Val")
        if idx<0: return ("NO_COL", h[:6])
        n=t=0
        for row in r:
            t+=1
            try:
                if float(row[idx])<0.05: n+=1
            except: pass
        return (n,t)

def nrows(path):
    with open(path, newline='') as f:
        return sum(1 for _ in csv.reader(f))-1

def count_true(path, col):
    with open(path, newline='') as f:
        r=csv.reader(f); h=next(r); idx=cidx(h,col)
        if idx<0: return ("NO_COL", h[:6])
        return sum(1 for row in r if len(row)>idx and row[idx].strip().upper()=="TRUE")

print("== DMP sig (adj.P.Val<0.05) ==")
for f in ["gse231574_DMP_Bor_vs_Ben.csv","gse231574_DMP_Mal_vs_Ben.csv","gse231574_DMP_Mal_vs_Bor.csv","gse231574_DMP_Trend.csv","gse179458_DMP_Mal_vs_Ben.csv"]:
    print(" ", f, sig_adjP(os.path.join(TAB,f)))
print("== DE sig (adj.P.Val<0.05) ==")
for f in ["gse78071_DE_PT_vs_FA.csv","gse78071_DE_Mal_vs_Ben.csv","gse78071_DE_Trend.csv","gse235267_DE_PT_vs_FA.csv","gse235267_DE_Mal_vs_Ben.csv"]:
    print(" ", f, sig_adjP(os.path.join(TAB,f)))
print("== network / GO ==")
print("  network edges:", nrows(os.path.join(TAB,"gse235267_lncRNA_mRNA_network.csv")))
print("  GO terms:", nrows(os.path.join(TAB,"gse235267_lncRNA_partner_GO.csv")))
print("== concordance rows ==", nrows(os.path.join(TAB,"gse78071_meth_module_concordance.csv")))
print("== druggable TRUE ==", count_true(os.path.join(TAB,"gse_pt_druggable_candidates.csv"),"druggable"))
print("== signature file ==")
sig = os.path.join(TAB,"gse78071_malignant_signature_score.csv")
if os.path.exists(sig):
    with open(sig,newline='') as f:
        for line in f: print("  "+line.rstrip()[:160])
else:
    print("  (skipped: gse78071_malignant_signature_score.csv no longer exists)")

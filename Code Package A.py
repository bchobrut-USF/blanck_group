from Bio.Seq import Seq
from Bio import pairwise2
import Bio.SeqIO as SeqIO
import pandas as pd
import regex as re
import numpy as np

def nt_matches(stringa, stringb):#returns tuple: (% match, nt/total string, True/False VJ > 20)
    u = zip(stringa, stringb)
    y=[]
    for i,j in u:
        if i==j:
            pass
        else: 
            y.append(j)
    try:
        percentj = float(len(stringa)-len(y))/float(len(stringb))
        percentj = percentj*100
        ntmatch = str((len(stringa)-len(y))) + "/" + str((len(stringb))) + " nt"
        vj = float(len(stringa)-len(y))
        return (percentj, ntmatch, vj)
    except:
        return (np.nan, np.nan, np.nan)
        
class VDJ(object):
    
    def score_only_align(self, subtype):
        if subtype == "V":
            ref_seqs = self.v_seqs
        elif subtype == "J":
            ref_seqs = self.j_seqs
        id_list = []
        alignments_list = []
        for ref in ref_seqs:
            alignments = pairwise2.align.localms(ref.seq, self.test_seq, 5, -10, -10, -10, score_only=True)
            alignments_list.append(alignments)
            id_list.append(ref.id)
        tmp = pd.DataFrame()
        tmp['%sID'%subtype] = id_list
        tmp['%sSCORE'%subtype] = alignments_list
        tmp = tmp.sort_values('%sSCORE'%subtype, ascending=False)
        tmp = tmp.reset_index(drop=True)
        tmp = tmp.iloc[0]
        return tmp
          
    def full_align(self, subtype):
        if subtype == "V":
            ref_seqs = self.v_seqs
        elif subtype == "J":
            ref_seqs = self.j_seqs
        id_list = []
        alignments_list = []
        for ref in ref_seqs:
            alignments = pairwise2.align.localms(ref.seq, self.test_seq, 5, -10, -10, -10, score_only=False)
            for alignment in alignments:
                alignments_list.append(alignment)
                id_list.append(ref.id)
        tmp = pd.DataFrame()
        tmp['%sID'%subtype] = id_list
        tmp['%sALIGNMENT'%subtype] = alignments_list
        tmp[['%sREFSEQ'%subtype, '%sTESTSEQ'%subtype, '%sSCORE'%subtype, '%sALIGNSTART'%subtype, '%sALIGNSTOP'%subtype]] = tmp['%sALIGNMENT'%subtype].apply(pd.Series)
        test_hit = lambda x : x['%sTESTSEQ'%subtype][x['%sALIGNSTART'%subtype]:x['%sALIGNSTOP'%subtype]]
        ref_hit = lambda x : x['%sREFSEQ'%subtype][x['%sALIGNSTART'%subtype]:x['%sALIGNSTOP'%subtype]]
        remainder_hit = lambda x : (x['%sTESTSEQ'%subtype][:x['%sALIGNSTOP'%subtype]])+(x['%sREFSEQ'%subtype][x['%sALIGNSTOP'%subtype]:])
        tmp['%sTESTALIGNSEQ'%subtype] = tmp.apply(test_hit, axis = 1)
        tmp['%sREFALIGNSEQ'%subtype] = tmp.apply(ref_hit, axis = 1)
        if subtype == "J":
            tmp['%sREMAINDERSEQ'%subtype] = tmp.apply(remainder_hit, axis = 1)
        tmp = tmp.sort_values('%sSCORE'%subtype, ascending=False)
        tmp = tmp.reset_index(drop=True)
        tmp = tmp.iloc[0]
        tmp["%s Match Percent"%subtype], tmp["%s Match"%subtype], tmp["%s Match Length"%subtype] = nt_matches(tmp['%sREFALIGNSEQ'%subtype], tmp['%sTESTALIGNSEQ'%subtype]) 
        return tmp
      
    def score_test(self):
        jdf = self.score_only_align("J")
        if jdf["JSCORE"] <= self.jthreshold:
            jdf['JID'] = "No Match"
            self.results = jdf
            self.results["JUNCTION"], self.results["CDR3"], self.results["REASON"] = "No Match", "No Match", "No Match"
            return False
        elif jdf["JSCORE"] >= self.jthreshold:
            vdf = self.score_only_align("V")
            if vdf["VSCORE"] <= self.vthreshold:
                vdf['VID'] = "No Match"
                self.results = pd.concat([vdf,jdf])
                self.results["JUNCTION"], self.results["CDR3"], self.results["REASON"] = "No Match", "No Match", "No Match"
                return False
            elif vdf["VSCORE"] >= self.vthreshold:
                return True
    
    def junction_find(self):
        junctions = re.search(r'((tgt|tgc)(...){4,20}(ttt|ttc|tgg)(ggt|ggc|gga|ggg)(...)(ggt|ggc|gga|ggg))', self.results['JREMAINDERSEQ'])
        if junctions is None:
            self.results["JUNCTION"], self.results["CDR3"], self.results["REASON"] = "Unproductive", "Unproductive", "Frame"
        else:
            junction = junctions.group(0)[:len(junctions.group(0))-9]
            translation = str(Seq(junction).ungap('-').translate())
            if '*' not in translation:
                self.results["JUNCTION"], self.results["CDR3"], self.results["REASON"] = junction, translation, "NA"
            else:
                self.results["JUNCTION"], self.results["CDR3"], self.results["REASON"] = "Unproductive", "Unproductive", "Stop"
 
 
    def run(self):
        score_test_pass = self.score_test()
        if score_test_pass == False:
            return self.results
        elif score_test_pass == True:
            vdf = self.full_align("V")
            jdf = self.full_align("J")
            self.results = pd.concat([vdf, jdf])
            self.junction_find()
            return self.results
 
 
class TRA(VDJ):
    def __init__(self, sequence, dbpath, vthreshold=65, jthreshold=65):
        self.test_seq = Seq(sequence).lower()
        self.v_seqs = list(SeqIO.parse("%sTRAV.fasta" %dbpath, "fasta"))
        self.j_seqs = list(SeqIO.parse("%sTRAJ.fasta" %dbpath, "fasta"))
        self.vthreshold = vthreshold
        self.jthreshold = jthreshold
        
class TRB(VDJ):
    def __init__(self, sequence, dbpath, vthreshold=65, jthreshold=65):
        self.test_seq = Seq(sequence).lower()
        self.v_seqs = list(SeqIO.parse("%sTRBV.fasta" %dbpath, "fasta"))
        self.j_seqs = list(SeqIO.parse("%sTRBJ.fasta" %dbpath, "fasta"))
        self.vthreshold = vthreshold
        self.jthreshold = jthreshold

class TRG(VDJ):
   def __init__(self, sequence, dbpath, vthreshold=65, jthreshold=65):
        self.test_seq = Seq(sequence).lower()
        self.test_seq = self.test_seq.reverse_complement()
        self.v_seqs = list(SeqIO.parse("%sTRGV.fasta" %dbpath, "fasta"))
        self.j_seqs = list(SeqIO.parse("%sTRGJ.fasta" %dbpath, "fasta"))
        self.vthreshold = vthreshold
        self.jthreshold = jthreshold

class TRD(VDJ):
    def __init__(self, sequence, dbpath, vthreshold=65, jthreshold=65):
        self.test_seq = Seq(sequence).lower()
        self.v_seqs = list(SeqIO.parse("%sTRDV.fasta" %dbpath, "fasta"))
        self.j_seqs = list(SeqIO.parse("%sTRDJ.fasta" %dbpath, "fasta"))
        self.vthreshold = vthreshold
        self.jthreshold = jthreshold
        
class IGH(VDJ):
    def __init__(self, sequence, dbpath, vthreshold=65, jthreshold=65):
        self.test_seq = Seq(sequence).lower()
        self.test_seq = self.test_seq.reverse_complement()
        self.v_seqs = list(SeqIO.parse("%sIGHV.fasta" %dbpath, "fasta"))
        self.j_seqs = list(SeqIO.parse("%sIGHJ.fasta" %dbpath, "fasta"))
        self.vthreshold = vthreshold
        self.jthreshold = jthreshold
    
class IGL(VDJ):
    def __init__(self, sequence, dbpath, vthreshold=65, jthreshold=65):
        self.test_seq = Seq(sequence).lower()
        self.v_seqs = list(SeqIO.parse("%sIGLV.fasta" %dbpath, "fasta"))
        self.j_seqs = list(SeqIO.parse("%sIGLJ.fasta" %dbpath, "fasta"))
        self.vthreshold = vthreshold
        self.jthreshold = jthreshold

class IGK(VDJ):
    def __init__(self, sequence, dbpath, vthreshold=65, jthreshold=65):
        self.test_seq = Seq(sequence).lower()
        self.v_seqs = list(SeqIO.parse("%sIGKV.fasta" %dbpath, "fasta"))
        self.j_seqs = list(SeqIO.parse("%sIGKJ.fasta" %dbpath, "fasta"))
        self.vthreshold = vthreshold
        self.jthreshold = jthreshold


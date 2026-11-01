"""Rebuild independent goldens: PYTHONPATH=<Biopython 1.86 install> python tools/reference_fixtures.py.
Biopython is a development reference only, not an application dependency.
"""
import itertools, json, random
from pathlib import Path
import Bio
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqUtils import MeltingTemp

root = Path(__file__).resolve().parents[1]
codons = [''.join(c) for c in itertools.product('ATCG', repeat=3)]
tables = {}
for table_id in (1, 2):
    table = CodonTable.unambiguous_dna_by_id[table_id]
    tables[str(table_id)] = {'translation': {c: str(Seq(c).translate(table=table_id)) for c in codons}, 'starts': table.start_codons, 'stops': table.stop_codons}
rng = random.Random(20261101)
sequences = ['ATGGCCATTGTAATGGGCCG','GTCGACCTGCAGCCCGGGAT','CGCGAATTCGCGATGGCCATTG','CGCGAATTCGCG','GCGCGCGC','ATGCATGCAT']
sequences += [''.join(rng.choice('ATCG') for _ in range(n)) for n in (15,18,20,22,24,26,28,30)]
tm = []
for sequence in sequences:
    selfcomp = sequence == str(Seq(sequence).reverse_complement())
    value = MeltingTemp.Tm_NN(sequence, nn_table=MeltingTemp.DNA_NN3, Na=50, saltcorr=5, dnac1=250 if selfcomp else 125, dnac2=0 if selfcomp else 125, selfcomp=selfcomp)
    tm.append({'sequence':sequence,'selfComplementary':selfcomp,'temperature':value})
result = {'reference':'Biopython','version':Bio.__version__,'seed':20261101,'conditions':{'Na_mM':50,'totalOligo_nM':250,'saltcorr':5,'nnTable':'DNA_NN3'},'sources':['https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi','https://github.com/biopython/biopython/blob/biopython-186/Bio/SeqUtils/MeltingTemp.py'],'tables':tables,'meltingTemperatures':tm}
out=root/'tests/fixtures/biopython-1.86.json';out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print('Wrote',out,'with',len(tm),'independent Tm cases and',len(codons)*2,'codon cases.')

#!/usr/bin/env python3
"""Validate scientific inputs, fitted results, figure assets, and portability."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
import pandas as pd
import pypdfium2 as pdfium
from docx import Document
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/reproduced';TARGET=ROOT/'data/validation_targets';DEST=ROOT/'validation';DEST.mkdir(exist_ok=True)
checks=[]
def check(name,condition,detail=''):
 checks.append(dict(check=name,status='PASS' if bool(condition) else 'FAIL',detail=str(detail)))
def compare(name,actual,expected,keys,values,tolerance=1e-6):
 a=pd.read_csv(actual);b=pd.read_csv(expected)
 x=a.merge(b,on=keys,suffixes=('_actual','_expected'),validate='one_to_one')
 check(name+' row identity',len(x)==len(a)==len(b),f'{len(a)} / {len(b)} rows')
 for v in values:
  va=x[v+'_actual'].to_numpy(float);vb=x[v+'_expected'].to_numpy(float)
  check(name+' '+v,np.allclose(va,vb,rtol=tolerance,atol=tolerance,equal_nan=True),float(np.nanmax(abs(va-vb))) if np.isfinite(va-vb).any() else 'NA')
 return a
mapping=pd.read_csv(ROOT/'data/functional_metagenomics/frozen_332_KO_functional_group_mapping.csv')
universe=pd.read_csv(ROOT/'data/functional_metagenomics/full_7931_KO_universe_review.csv')
check('Unique primary KO assignments',len(mapping)==332 and mapping.ko_id.nunique()==332 and mapping.functional_group.nunique()==13)
check('Complete KO universe',len(universe)==7931 and universe.ko_id.nunique()==7931 and set(mapping.ko_id)<=set(universe.ko_id))
amp=pd.read_csv(OUT/'amplicons/amplicon_global_statistics.csv')
expected={('field','ITS'):(18,416,416),('field','16S'):(19,961,960),('microcosm','ITS'):(40,240,240),('microcosm','16S'):(36,634,634)}
for (experiment,marker),(N,pre,post) in expected.items():
 a=amp[(amp.experiment==experiment)&(amp.marker==marker)].iloc[0]
 check(f'{experiment} {marker} sample and ASV counts',(int(a.N),int(a.taxa_after_prevalence),int(a.taxa_after_rarefaction))==(N,pre,post))
 check(f'{experiment} {marker} prevalence floor',int(a.minimum_prevalence_samples)==max(2,int(np.ceil(.05*N))) and a.minimum_observed_asv_prevalence>=2)
 if experiment=='field':
  for suffix in ['filtered_ASV_counts','ASV_relative_abundance','genus_relative_abundance','genus_composition_summary','phylum_composition_summary']:
   g=pd.read_csv(OUT/f'amplicons/field_{marker}_{suffix}.csv');f=pd.read_csv(ROOT/f'data/amplicons/derived/field_{marker}_{suffix}.csv')
   check(f'Field {marker} derived {suffix}',g.shape==f.shape and g.equals(f))
 else:
  a=pd.read_csv(OUT/f'amplicons/microcosm/{marker}/sequential_permanova.csv');b=pd.read_csv(TARGET/'microcosm_amplicon_permanova.csv');b=b[b.marker==marker]
  for c in ['F','p','R2']:check(marker+' factorial '+c,np.allclose(a[c],b[c],rtol=1e-7,atol=1e-9,equal_nan=True))
  disp=pd.read_csv(OUT/f'amplicons/microcosm/{marker}/dispersion_tests.csv');ref=pd.read_csv(TARGET/'microcosm_amplicon_dispersion.csv');ref=ref[ref.dataset.eq('Microcosm '+marker)]
  check(marker+' PERMDISP',np.allclose(disp.p,ref.p,atol=1e-12) and np.allclose(disp.F,ref.F,atol=1e-9))
  partition=pd.read_csv(OUT/f'amplicons/microcosm/{marker}/sequential_variance_partition.csv')
  check(marker+' variance sum',abs(partition.proportion.sum()-1)<1e-12)
  order=pd.read_csv(OUT/f'amplicons/microcosm_{marker}_composition_sample_order.csv')
  origins=['InocS+','InocUV-','InocS-'];substrates=['SusS+','SusUV-','SusS-']
  tuples=[(substrates.index(r.Treatment),origins.index(r.Inoculum2),int(r.Block),r.sample_id) for r in order.itertuples()]
  check(marker+' composition order',tuples==sorted(tuples) and len(order)==N)
field=amp[amp.experiment=='field'];ref=pd.read_csv(TARGET/'field_amplicon_statistics.csv')
for c in ['F','R2','p','dispersion_p']:
 z=field.merge(ref,on='marker',suffixes=('_a','_b'));check('Field amplicon '+c,np.allclose(z[c+'_a'],z[c+'_b'],atol=1e-8))
compare('Field functional profile',OUT/'functional_metagenomics/field/global_statistics.csv',TARGET/'functional_field_statistics.csv',['experiment'],['N','F','R2','p','dispersion_p'])
compare('Microcosm functional profile',OUT/'functional_metagenomics/microcosm/sequential_permanova.csv',TARGET/'functional_microcosm_permanova.csv',['term'],['F','R2','p'])
groups=compare('Functional-group models',OUT/'functional_metagenomics/microcosm/factorial_anova.csv',TARGET/'functional_group_anova.csv',['functional_group','term'],['F','p','q_BH'])
measure=compare('Measurement models',OUT/'functional_measurements/factorial_anova.csv',TARGET/'measurements_anova.csv',['outcome','term'],['N','numerator_df','denominator_df','F','p','p_BH'])
compare('Measurement substrate contrasts',OUT/'functional_measurements/substrate_pairwise_Tukey.csv',TARGET/'measurement_substrate_contrasts.csv',['outcome','group1','group2'],['estimate','p_Tukey'])
measurement_input=pd.read_csv(ROOT/'data/functional_measurements/cleaned_functional_measurements.csv')
inoculated=measurement_input[~measurement_input.is_control].set_index('id')
expected_bottles=set(range(1,46))-{29,35}
check('Measurement input bottle identities',measurement_input.id.is_unique and set(inoculated.index)==expected_bottles)
check('Measurement input excludes block 6',not measurement_input.block.eq(6).any() and set(inoculated.block)==set(range(1,6)))
check('Measurement input controls retained',set(measurement_input.loc[measurement_input.is_control,'id'])==set(range(49,56)))
check('Measurement inclusion ledger',set(pd.read_csv(ROOT/'data/sample_inclusion_ledger.csv').query("assay == 'functional_measurements' and included_primary").sample_id.astype(int))==expected_bottles)
plotted=pd.read_csv(OUT/'functional_measurements/plotted_samples.csv')
measurement_columns={'decomposition':'decomposition_pct','biomass':'microbial_biomass_c_ug_g','beta_glucosidase':'beta_glucosidase','cue_proxy':'cue_proxy'}
for outcome,column in measurement_columns.items():
 expected_ids=expected_bottles-({9,44,45} if outcome=='cue_proxy' else set())
 points=plotted[plotted.outcome==outcome].set_index('sample_id')
 check(outcome+' available measurement IDs',set(inoculated.index[inoculated[column].notna()])==expected_ids)
 check(outcome+' plotted bottle IDs',points.index.is_unique and set(points.index)==expected_ids)
 check(outcome+' plotted values',np.allclose(points.value,inoculated.loc[points.index,column],rtol=1e-12,atol=1e-12))
 check(outcome+' fitted sample size',measure.loc[measure.outcome==outcome,'N'].eq(len(expected_ids)).all())
design=pd.read_csv(OUT/'functional_measurements/design_checks.csv')
check('Measurement random block design',len(design)==4 and design.blocks.eq(5).all() and design.fixed_columns.eq(9).all() and design.fixed_block_columns.eq(0).all() and not design.singular.any())
brackets=pd.read_csv(OUT/'functional_measurements/displayed_brackets.csv')
contrasts=pd.read_csv(OUT/'functional_measurements/substrate_pairwise_Tukey.csv')
significant=contrasts[contrasts.p_Tukey.lt(.05)]
check('Measurement displayed brackets',set(zip(brackets.outcome,brackets.group1,brackets.group2))==set(zip(significant.outcome,significant.group1,significant.group2)))
compare('Genus overlap',OUT/'amplicon_comparison/data/venn_counts.csv',TARGET/'genus_overlap_counts.csv',['marker'],['field_only','shared','microcosm_only'])
for frame,qcol,n,label in [(groups,'q_BH',12,'functional group'),(measure,'p_BH',4,'measurement')]:
 for term,g in frame.groupby('term',sort=False):
  values=g.p.to_numpy();idx=np.argsort(values);q=np.minimum.accumulate((values[idx]*n/np.arange(1,n+1))[::-1])[::-1].clip(0,1)
  check(label+' BH '+term,len(g)==n and np.allclose(g[qcol].to_numpy()[idx],q,atol=1e-12))
check('Four measurement substrate effects',measure[measure.term=='substrate'].p_BH.lt(.05).all())
check('Measurement inoculum and interaction tests',measure[measure.term!='substrate'].p_BH.ge(.05).all())
for table,n in [('Table_S1_functional_group_definitions.csv',13),('Table_S2_complete_KO_mapping.csv',332)]:check(table,len(pd.read_csv(ROOT/'tables'/table))==n)
assets=json.loads((ROOT/'docs/figure_manifest.json').read_text())
for a in assets:
 f=ROOT/'figures/final'/a['file'];check(a['file']+' checksum',f.exists() and hashlib.sha256(f.read_bytes()).hexdigest()==a['sha256'])
 if f.suffix=='.pdf':
  doc=pdfium.PdfDocument(f);check(a['file']+' single page',len(doc)==1);doc.close()
check('Final figure set',{x['file'] for x in assets}=={'Figure1.pptx',*[f'Figure{i}.pdf' for i in range(2,7)],*[f'FigS{i}.pdf' for i in range(1,12)]})
doc=Document(ROOT/'figures/final/Supplementary-Figures.docx')
check('Supplementary bundle has eleven figures',len(doc.inline_shapes)==11 and sum(p.text.startswith('Figure S') for p in doc.paragraphs)==11)
check('Supplementary table bundle',len(Document(ROOT/'tables/Supplementary_Tables_S1-S2.docx').tables)==2)
for folder in ['scripts','docs']:
 for p in (ROOT/folder).rglob('*'):
  if p.is_file() and p.suffix.lower() in ['.r','.py','.sh','.md','.json','.txt']:
   s=p.read_text();check('Portable '+str(p.relative_to(ROOT)),('/home/'+'felius') not in s and ('wsl'+'.localhost') not in s)
failed=[c for c in checks if c['status']=='FAIL']
report={'passed':not failed,'checks':checks}
(DEST/'VALIDATION_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
(DEST/'VALIDATION_REPORT.txt').write_text('\n'.join(f"[{c['status']}] {c['check']}: {c['detail']}" for c in checks)+'\n')
manifest=[]
for p in sorted(ROOT.rglob('*')):
 if p.is_file() and not any(x in p.relative_to(ROOT).parts for x in ['.git','results','validation','__pycache__']):manifest.append(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT)}')
(DEST/'MANIFEST.sha256').write_text('\n'.join(manifest)+'\n')
print(f'{len(checks)} release checks; {len(failed)} failed.')
if failed:
 for c in failed:print(c)
 raise SystemExit(1)
print('All release checks passed.')

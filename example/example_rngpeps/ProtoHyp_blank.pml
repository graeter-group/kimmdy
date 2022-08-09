load pep_i.pdb, pep_i
edit resn PRO and resid pos and name 2HG
replace O,4,2
alter pep_i and resid pos,resn = "HYP"
alter pep_i and resid pos and name CD, name = "CD2"
alter pep_i and resid pos and name O01, name = "OD1"
set pdb_use_ter_records, 0
save pep_i.pdb, pep_i

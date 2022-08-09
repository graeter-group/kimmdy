load pep_i.pdb, pep_i
edit resn TYR and resid pos and name 1HE
replace O,4,2
alter pep_i and resid pos,resn = "DOP"
alter pep_i and resid pos and name OH, name = "OH1"
alter pep_i and resid pos and name O01, name = "OH2"
set pdb_use_ter_records, 0
save pep_i.pdb, pep_i

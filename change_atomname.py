# Script to modify .gro and .top in order to use dissociation energies specific in crosslinks (L4Y and L5Y) or glycines
# creates new atom names for these residues where we calculated deviating dissociation energies


import functions as func

topfile = 'topol.top'
newtop_file = 'topol_modified.top'
grofile = 'lin_comm_grps2018.gro'
newgro_file = 'gro_modified.gro'

top_all, top_array = func.get_data_from_file(topfile)
gro_all, gro_array = func.get_data_from_file(grofile)

residues = ['L4Y', 'L5Y']  # where to make changes


for i in range(len(top_all)):

     if 'GLY' in top_all[i]:
        if 'CA' in top_all[i]:
            top_all[i] = top_all[i].replace(' CA', 'CAG')

    if 'L4Y' in top_all[i]:
        if 'CD' in top_all[i]:
            top_all[i] = top_all[i].replace(' CD', 'CD4')
        if 'CE' in top_all[i]:
            top_all[i] = top_all[i].replace(' CE', 'CE4')

    if 'L5Y' in top_all[i]:
        if 'CD' in top_all[i]:
            top_all[i] = top_all[i].replace(' CD', 'CD5')
        if 'CE' in top_all[i]:
            top_all[i] = top_all[i].replace(' CE', 'CE5')
            
func.store_linelist_to_file(top_all, newtop_file)

for i in range(len(gro_all)):

     if 'GLY' in gro_all[i]:
        if 'CD' in gro_all[i]:
            gro_all[i] = gro_all[i].replace(' CA', 'CAG')

     if 'L4Y' in gro_all[i]:
        if 'CD' in gro_all[i]:
            gro_all[i] = gro_all[i].replace(' CD', 'CD4')
        if 'CE' in gro_all[i]:
            gro_all[i] = gro_all[i].replace(' CE', 'CE4')

     if 'L5Y' in gro_all[i]:
        if 'CD' in gro_all[i]:
            gro_all[i] = gro_all[i].replace(' CD', 'CD5')
        if 'CE' in gro_all[i]:
            gro_all[i] = gro_all[i].replace(' CE', 'CE5')

func.store_linelist_to_file(gro_all, newgro_file)

        

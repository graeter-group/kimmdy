#%%
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from collections.abc import Iterable
from typing import Generator
from copy import deepcopy
import itertools
from operator import itemgetter

#%%
# from parsing
Topology = dict[str, list[list[str]]]

def get_sections(
    seq: Iterable[str], section_marker: str
) -> Generator[list[str], None, None]:
    data = [""]
    for line in seq:
        if line.startswith(section_marker):
            if data:
                # first element will be empty
                # because newlines mark sections
                data.pop(0)
                # only yield section if non-empty
                if data:
                    yield data
                data = []
        data.append(line.strip("\n"))
    if data:
        yield data

def extract_section_name(ls: list[str]) -> tuple[str, list[str]]:
    """takes a list of lines and return a tuple
    with the name and the lines minus the
    line that contained the name.
    Returns the empty string of no name was found.
    """
    for i, l in enumerate(ls):
        if l and l[0] != ";" and "[" in l:
            name = l.strip("[] \n")
            ls.pop(i)
            return (name, ls)
    else:
        return ("", ls)

def read_topol(path: Path) -> Topology:
    # TODO look into following #includes
    # TODO look into [ intermolecule ] section
    with open(path, "r") as f:
        sections = get_sections(f, "\n")
        d = {}
        for i, s in enumerate(sections):
            # skip empty sections
            if s == [""]:
                continue
            name, content = extract_section_name(s)
            content = [c.split() for c in content if c]
            if not name:
                name = f"BLOCK {i}"
            # sections can be duplicated.
            # append values in such cases
            if name not in d:
                d[name] = content
            else:
                d[name] += content
        return d

def write_topol(d: Topology, outfile: Path) -> None:
    with open(outfile, "w") as f:
        for title, content in d.items():
            if title.startswith("BLOCK "):
                f.write(f"\n")
            else:
                f.write(f"[ {title} ]\n")
            s = (
                "\n".join(
                    [
                        " ".join([x.ljust(8) if l[0] != ";" else x for x in l])
                        for l in content
                    ]
                )
                + "\n\n"
            )
            f.write(s)

class ConversionType(Enum):
    BREAK = auto()
    MOVE = auto()

@dataclass
class ConversionRecipe:
    """A ConversionRecipe.
    encompasses a single transformation, e.g. moving one
    atom or braking one bond.
    Parameters
    ----------
    type : list[ConversionType.BREAK or .MOVE]
    atom_idx : list[(from, to)]
    """

    type: list[ConversionType] = field(default_factory=list)
    atom_idx: list[tuple[int, int]] = field(default_factory=list)

#%%
# new implementation


class Atom:
    def __init__(self,idx: str,atomtype = None, atomname = None, resname = None):
        self.idx = idx
        self.atomtype = atomtype
        self.atomname = atomname
        self.resname = resname
        self.bound_to = []

def str_to_int(elem):
    try:
        return(int(elem))
    except ValueError:
        raise ValueError("Not all XList elements are integers!")

def eval_bond(entry):
    return((str_to_int(entry[0]),str_to_int(entry[1])))

def eval_angle(entry):
    return((str_to_int(entry[1]),str_to_int(entry[0]),str_to_int(entry[2])))

def eval_dihedral(entry):
    return((str_to_int(entry[1]),str_to_int(entry[2]),str_to_int(entry[0]),str_to_int(entry[3])))

def check_idx(object):
    try:
        return(str_to_int(object.idx))
    except:
        raise ValueError("Non Atom object in AtomList")

class localGraph:
    def __init__(self,topoldict,heavy_idx,ffdir):
        self.topoldict = topoldict
        self.heavy_idx = heavy_idx
        self.radicals = []
        self.ffdir = ffdir
        self.ff = {}

        self.AtomList = [Atom(self.heavy_idx)]
        self.BondList = []
        self.PairList = []
        self.AngleList = []
        self.DihedralList = []

    def order_lists(self):
        self.AtomList = sorted(self.AtomList,key=check_idx)   #doesn't work with Atom class anymore
        self.BondList = sorted(self.BondList,key=eval_bond)
        self.PairList = sorted(self.PairList,key=eval_bond)
        self.AngleList = sorted(self.AngleList,key=eval_angle)
        self.DihedralList = sorted(self.DihedralList,key=eval_dihedral)

    def get_AtomList_property(self,property):
        if property == 'idx':
            return [x.idx for x in self.AtomList]
        elif property == 'bound_to':
            return [x.bound_to for x in self.AtomList]
        elif property == 'atomtype':
            return [x.atomtype for x in self.AtomList]
        elif property == 'atomname':
            return [x.atomname for x in self.AtomList]
        elif property == 'resname':
            return [x.resname for x in self.AtomList]
        else:
            raise ValueError("Can't find property {} in object".format(property))  

    def update_bound_to(self):
        atoms_idx = self.get_AtomList_property('idx')
        bound_dict = {}
        for atom_idx in atoms_idx:
            bound_dict[atom_idx] = []
        for bond in self.BondList:
            if not bond[1] in bound_dict[bond[0]]:
                bound_dict[bond[0]].append(bond[1])
            if not bond[0] in bound_dict[bond[1]]:
                bound_dict[bond[1]].append(bond[0])
        for atom_idx in atoms_idx:
            self.AtomList[atoms_idx.index(atom_idx)].bound_to = bound_dict[atom_idx]

    def construct_graph(self, depth = 3):
        """
        searches the bonds section of a topology up to depth bonds deep
        to create a local graph i.e. filling the AtomList and BondList
        """
        curratoms = [self.heavy_idx]
        for i in range(depth):
            addlist = []
            for bond in self.topoldict['bonds']:
                if any(x in bond[:2] for x in curratoms):
                    if not bond[:2] in self.BondList:
                        self.BondList.append(bond[:2])
                    for j in [0,1]:
                        if not bond[j] in [x.idx for x in self.AtomList]:
                            addlist.append(bond[j])
                            self.AtomList.append(Atom(bond[j]))
            #print(f"{i}, current atoms: {curratoms}, addlist: {addlist}, AtomList: {self.AtomList}")
            curratoms = deepcopy(addlist)
        
        self.order_lists()
        self.update_bound_to()

    def build_PADs(self):
        """
        build pairs, angles dihedrals from BondList
        in an ordered fashion
        """
        atom_idx = self.get_AtomList_property('idx')
        # defining angles
        for atom in self.AtomList:
            if len(atom.bound_to) >= 2:
                print(atom.idx,itertools.combinations(atom.bound_to,2))
                for comb in itertools.combinations(atom.bound_to,2):
                    self.AngleList.append([comb[0],atom.idx,comb[1]])

        # defining dihedrals
        for atom1 in self.AtomList:
            if len(atom1.bound_to) >= 2:
                for partner in atom1.bound_to:
                    atom2 = self.AtomList[atom_idx.index(partner)]
                    if len(atom2.bound_to) >= 2 and int(atom2.idx) > int(atom1.idx):
                        atom1_periphery = [x for x in atom1.bound_to if x != atom2.idx]
                        atom2_periphery = [x for x in atom2.bound_to if x != atom1.idx]
                        for prod in itertools.product(atom1_periphery,atom2_periphery):
                            self.DihedralList.append([prod[0],atom1.idx,atom2.idx,prod[1]])

        # all pairs have a dihedral between them
        for dihedral in self.DihedralList:
            if int(dihedral[0]) < int(dihedral[3]):
                self.PairList.append([dihedral[0],dihedral[3]])
            else:
                self.PairList.append([dihedral[3],dihedral[0]])

        self.order_lists()

    def search_terms(self,TypeList,idx:str,center):
        CopyList = deepcopy(TypeList)
        hits = []
        if center and len(TypeList[0]) > 2:
            start = 1
            end = -1
        else:
            start = 0
            end = None
        for entry in CopyList:
            if idx in entry[start:end]: 
                hits.append(entry)
        return hits

    def get_terms_with_atom(self,idx:str,center = False, add_function = False):
        #print(self.get_AtomList_property('idx'))
        if not idx in self.get_AtomList_property('idx'):
            return

        termdict = {}
        termdict["bonds"] = self.search_terms(self.BondList,idx,center)
        termdict["pairs"] = self.search_terms(self.PairList,idx,center)
        termdict["angles"] = self.search_terms(self.AngleList,idx,center)
        termdict["dihedrals"] = self.search_terms(self.DihedralList,idx,center)

        if add_function:
            termdict_funct = deepcopy(termdict)
            funct_dict = {'bonds':'1','pairs':'1','angles':'1','dihedrals':'9'}
            for section in ['bonds','pairs','angles','dihedrals']:
                termdict_funct[section] = [[*x,funct_dict[section]] for x in termdict_funct[section]]
            return termdict_funct
        return termdict

    def add_Atom(self,atom: Atom):
        if not atom in self.AtomList:
            self.AtomList.append(atom)
        else: print('Atom already exists!')

    def add_Bond(self,bond: list):
        if all(x in self.get_AtomList_property('idx') for x in bond):
            self.BondList.append(bond)
            self.update_bound_to()
        else: print('Could not establish bond!')

    def remove_terms(self,termdict):
        graphdict = {'bonds':self.BondList,'pairs':self.PairList,'angles':self.AngleList,'dihedrals':self.DihedralList}
        for section in ['pairs','bonds','angles','dihedrals']:
            try:
                for term in termdict[section]:
                    graphdict[section].remove(term)
                    print(f"removed term {term} from graph {section}")
            except ValueError:
                print(f"Couldn't find term {term} in TypeList")

        self.update_bound_to()

    def compare_section(self,section,TypeList):
        """
        checks whether all entries in TypeList are contained in the 
        respective section.
        """
        typelen = len(TypeList[0])
        missing_terms = deepcopy(TypeList)
        for entry in section:
            #print(missing_terms)
            if entry[:typelen] in missing_terms:
                missing_terms.remove(entry[:typelen])
        return missing_terms

    def find_missing_terms(self):
        termdict = {}
        termdict["bonds"] = self.compare_section(self.topoldict['bonds'],self.BondList)
        termdict["pairs"] = self.compare_section(self.topoldict['pairs'],self.PairList)
        termdict["angles"] = self.compare_section(self.topoldict['angles'],self.AngleList)
        termdict["dihedrals"] = self.compare_section(self.topoldict['dihedrals'],self.DihedralList)
        return termdict

    def get_atomprop(self,property):
        atom_idx = self.get_AtomList_property('idx')
        for entry in self.topoldict['atoms']:
            if entry[0] in atom_idx:
                if property == 'atomtype':
                    self.AtomList[atom_idx.index(entry[0])].atomtype = entry[1]
                elif property == 'atomname':
                    self.AtomList[atom_idx.index(entry[0])].atomname = entry[4]
                elif property == 'resname':
                    self.AtomList[atom_idx.index(entry[0])].resname = entry[3]
                else:
                    return

    def is_radical(self,atom_idx):
        print('\n\n\n')
        if None in self.get_AtomList_property('atomtype'):
            self.get_atomprop('atomtype')

        bonds = {'H':1,'HC':1,'H1':1,'O':1,'N':3,'C':3,'CT':4} #not exhaustive,yet 

        atoms_idx = self.get_AtomList_property('idx')
        atoms_atomtype = self.get_AtomList_property('atomtype')    
        print(self.AtomList[atoms_idx.index(atom_idx)].bound_to)
        if len(self.AtomList[atoms_idx.index(atom_idx)].bound_to) < bonds[atoms_atomtype[atoms_idx.index(atom_idx)]]:
            print(f"{atom_idx} is a radical")
            return True
        else:
            print(f"{atom_idx} is not a radical")
            return False


    def find_radical(self):
        if None in self.get_AtomList_property('atomtype'):
            self.get_atomprop('atomtype')

        bonds = {'H':1,'HC':1,'H1':1,'O':1,'N':3,'C':3,'CT':4} #not exhaustive,yet 
        radicals = []

        atom_idx = self.get_AtomList_property('idx')
        atoms_atomtype = self.get_AtomList_property('atomtype')    
        for i in range(1,len(atom_idx[:-1])):            #[1:-1] because there are bonds missing at the ends of the localgraph
            if len(self.AtomList[i].bound_to) < bonds[atoms_atomtype[i]]:
                radicals.append(atom_idx[i])
        self.radicals = radicals

    def get_ff_sections(self):
        ffbonded = read_topol(self.ffdir / "ffbonded.itp")          #fringe use?!
        self.ff['bondtypes'] = ffbonded['bondtypes']
        self.ff['angletypes'] = ffbonded['angletypes']

    def parameterize_prop_terms(self,prop,terms):
        #print(terms)
        terms = deepcopy(terms)
        terms_prm = []

        atom_idx = self.get_AtomList_property('idx')
        atoms_atomtype = self.get_AtomList_property('atomtype')    

        terms_mapped = [[atoms_atomtype[atom_idx.index(x)] for x in term] for term in terms]
        #print(terms_mapped)
        
        for i,term in enumerate(terms_mapped):
            for entry in self.ff[prop]:
                if term == entry[:len(term)] or term[::-1] == entry[:len(term)]:
                    stop = entry.index(';')
                    #print(entry[len(term):stop])
                    terms_prm.append([*terms[i],*entry[len(term):stop]])
                    break
            else:
                print(f"\n\n No parameters found for {term}/{terms[i]}!\n\n")
        #print(terms_prm)
        return terms_prm

    def patch_bond(self,bonds,atom_idx):
        newfrac = 0.98  # factor by which all non-hydrogen bonds of the atom_idx are corrected
        SOfrac = 0.955

        atoms_idx = self.get_AtomList_property('idx')
        atoms_atomtype = self.get_AtomList_property('atomtype')
        atoms_atomname = self.get_AtomList_property('atomname')

        for bond in bonds:                    
            partnerpos = 0 if str(atom_idx) == bond[1] else 1
            radicalpos = bond[1-partnerpos]
            partner_idx = bond[partnerpos]
            partner_atomtype = atoms_atomtype[atoms_idx.index(partner_idx)]
            radical_atomtype = atoms_atomtype[atoms_idx.index(atom_idx)]
            radicalname = atoms_atomname[atoms_idx.index(atom_idx)]
            
            #print(partner_idx,partner_atomtype,atom_idx,radical_atomtype,radicalname)

            if partner_atomtype.startswith('H'):
                #print(f"{partner_idx} is a hydrogen, no change of parameters necessary!")
                continue

            req = float(bond[3]) * newfrac      # change to r_eq happens here
            k = float(bond[4]) 
            #print(req,k)
            
            #special fixes to bond length that are added on top of the general *0.98 factor
            if partner_atomtype == 'N' and radical_atomtype == 'CT' and radicalname == 'CA':      # N-CA of C-alpha atom_idx
                print("!!Found CaR!!")
                req -= 0.006

            if partner_atomtype == 'C' and radical_atomtype == 'N' and radicalname == 'N':        # C-N of N atom_idx
                print("!!Found bb NR!!")
                req += 0.008

            if any(at in ['S','O'] for at in [radical_atomtype,partner_atomtype]):               
                req = req * (SOfrac/newfrac) # correcting by a different value

            bond[3] = '{:7.5f}'.format(req)         # will carry back to atom_terms 
            bond[4] = '{:13.6f}'.format(k)
            bond.append(' ; patched parameter')

    def patch_angle(self,angles,atom_idx):
        aromatic_offset = 10    # add this to all angle theteq with a center atom_idx aromatic atom type
        newtheteq = 117         # new theteq of all angles around the center atom_idx

        atoms_idx = self.get_AtomList_property('idx')
        atoms_atomtype = self.get_AtomList_property('atomtype')
        atoms_atomname = self.get_AtomList_property('atomname')

        for angle in angles:
            radical_atomtype = atoms_atomtype[atoms_idx.index(atom_idx)]

            theteq = float(angle[4])
            k = float(angle[5])
            #print(k,theteq)

            if radical_atomtype in ['CR','NA','CW','CC','CA','CN','C*']:
                theteq += aromatic_offset
            else:
                theteq = newtheteq

            angle[4] = '{:11.7f}'.format(theteq)        # will carry back to atom_terms 
            angle[5] = '{:10.6f}'.format(k)
            angle.append(' ; patched parameter')

    def patch_dihedral_CA(self,dihedrals,atom_idx):
        phivals = ['1.6279944','21.068532', '1.447664']         # from own MCSA
        psivals = ['6.556746', '20.284450', '0.297901']

        atoms_idx = self.get_AtomList_property('idx')
        atoms_atomtype = self.get_AtomList_property('atomtype')
        atoms_atomname = self.get_AtomList_property('atomname')
        atom_resname = self.get_AtomList_property('resname')

        radical_resname = atom_resname[atoms_idx.index(atom_idx)]
        radical_atomname = atoms_atomname[atoms_idx.index(atom_idx)]

        if radical_resname in ['Pro','Hyp'] or not radical_atomname == 'CA':
            print(f"Did not add a dihedral term because the atom_idx is in Pro or Hyp (here {radical_resname}) or not a C-alpha (here {radical_atomname})")
            return []

        phi = []
        psi = []
        dihedrals_mapped = [[atoms_atomname[atoms_idx.index(x)] for x in dihedral] for dihedral in dihedrals]
        for i,dihedral in enumerate(dihedrals_mapped):
            print(i,dihedral)
            if dihedral == ['C','N','CA','C']:      #Phi
                print('Phi: ',dihedrals[i])
                for ii in range(3):
                    phi.append([*dihedrals[i],'9','180.000000',phivals[ii],str(ii+1),' ; patched parameter'])

            elif dihedral == ['N','CA','C','N']:     #Psi
                print('Psi: ',dihedrals[i])
                for ii in range(3):
                    psi.append([*dihedrals[i],'9','180.000000',psivals[ii],str(ii+1),' ; patched parameter'])


        return [*phi,*psi]

    def patch_improper(self,dihedrals,atom_idx):
        newphik = '4.6024000'        # same as backbone N improper
        atoms_idx = self.get_AtomList_property('idx')
        radical_Atom = self.AtomList[atoms_idx.index(atom_idx)]

        print(f"atom_idx is bound to {len(radical_Atom.bound_to)} Atoms.")
        if len(radical_Atom.bound_to) == 3:         #this would mean atom_idx went from 4 to 3 partners -> tetrahedral to planar
            improper = [radical_Atom.bound_to[0],radical_Atom.bound_to[1],str(atom_idx),radical_Atom.bound_to[2],'4','180.0000000',newphik,'2',' ; patched parameter']     #improper entry
            return improper
        else:
            return None

    def parameterize_around_atom(self,atom_idx):
        if None in self.get_AtomList_property('atomtype'):
            self.get_atomprop('atomtype')
            self.get_atomprop('atomname')
            self.get_atomprop('resname')
        print('\n --- Parameterizing Radicals ---')
        print(self.get_AtomList_property('idx'))
        print(self.get_AtomList_property('atomtype'))
        print(self.get_AtomList_property('atomname'))
        print(self.get_AtomList_property('resname'))
        #self.find_radical()
        #print(self.radicals)
        print(self.AngleList)
  
        atom_terms = self.get_terms_with_atom(atom_idx,center=True)
        # don't need pairs (and dihedrals)
        atom_terms['pairs'].clear()
        #del atom_terms['dihedrals']
        # deal with ff
        print('a1',atom_terms)
        self.get_ff_sections()
        atom_terms['bonds'] = self.parameterize_prop_terms('bondtypes',atom_terms['bonds'])
        atom_terms['angles'] = self.parameterize_prop_terms('angletypes',atom_terms['angles'])
        print('a2',atom_terms)
        # apply patches
        if self.is_radical(atom_idx):
            self.patch_bond(atom_terms['bonds'],atom_idx)
            self.patch_angle(atom_terms['angles'],atom_idx)
            atom_terms['dihedrals'] = self.patch_dihedral_CA(atom_terms['dihedrals'],atom_idx)
            atom_terms['improper'] = self.patch_improper(atom_terms['dihedrals'],atom_idx)

        
        print('atom_idx terms',atom_terms)
        return atom_terms
            

def find_heavy(bonds,H_idx: str):
    """
    takes the bonds section of a topology object and an atom index
    to return the index of the first found bond partner.
    """
    for bond in bonds:
        if H_idx in bond[:2]:
            if H_idx == bond[0]:
                return bond[1]
            elif H_idx == bond[1]:
                return bond[0]
    raise ValueError('Atom from ConversionRecipe not found in "bonds" section of the topology')

def compare_ge(term1,term2,size):           #typically, term1 is from the topology, term2 will be put in the topology
    term1_int = []
    term2_int = []
    orderdict = {2:[0,1],3:[1,0,2],4:[1,2,0,3]}
    for i in range(size):
        term1_int.append(str_to_int(term1[i]))
        term2_int.append(str_to_int(term2[i]))

    if term1_int[orderdict[size][0]] > term2_int[orderdict[size][0]]:       #always insert term when main index is higher in the insert term
        print(term1_int[orderdict[size][0]], term2_int[orderdict[size][0]])
        print('True1')
        return True
    for i in orderdict[size]:
        #print(term1_int,term2_int)
        #print('compare',term1_int[i],term2_int[i])
        #print(term1_int[i] >= term2_int[i])
        if term1_int[i] > term2_int[i]:
            print('True2')
            return True                    # False means term2 will not be inserted before term1
        elif term1_int[i] == term2_int[i]:
            continue
    
        return False                              

def topol_add_terms(topoldict,adddict):
    print(f"\naddict: {adddict}")
    copydict = deepcopy(topoldict)
    insertdict = {'pairs':[],'bonds':[],'angles':[],'dihedrals':[]} 

    mapsection = {'pairs':2,'bonds':2,'angles':3,'dihedrals':4}            # number of elements we want to compare
    for section in ['pairs','bonds','angles','dihedrals']:     
        k = mapsection[section]
        start = 0
        for term in adddict[section]:
            
            for i,term_topol in enumerate(topoldict[section][start:],start):
                if len(term_topol) != 0 and term_topol[0] != ';':
                    if term[:k] == term_topol[:k]:
                        copydict[section][i] = term
                        start = i
                        print('added term: ', term)
                        break   
                    elif compare_ge(term_topol,term,k):             #not sure with this
                        print('ge TRUE',term_topol,term)
                        insertdict[section].append([term,i])
                        start = i
                        break
            else:
                print(f"\nDid not find a spot for term {term}!")
        print(f"section {section} inserts: {insertdict[section]}")
        for entry in insertdict[section][::-1]:
            print('entry:',entry)
            copydict[section].insert(entry[1],entry[0])

        # get rid of explicit dihedrals after implicit dihedrals / this step may be unnecessary if we find a better way to go from explicit dihedral to implicit dihedral
        rmv_dih = []
        prevdih = []
        #print('copydict',copydict['dihedrals'])
        for i,entry in enumerate(copydict['dihedrals']):
            if entry[:4] == prevdih[:4]:
                #print(len(prevdih),len(entry))
                if len(prevdih) == 5 and len(entry) > 5:                   
                    rmv_dih.append(i)
                    continue
            else:
                prevdih = entry

    print('removing dihedrals at these positions:', rmv_dih)
    for val in rmv_dih[::-1]:
        del copydict['dihedrals'][val]

    if adddict['improper']:                                 #ugly construct, maybe reverse search because impropers are at the end of dihedral
        for i,term_topol in enumerate(topoldict['dihedrals']):
            if adddict['improper'][:5] == term_topol[:5]:       
                print('added term: ', adddict['improper'])     
                copydict['dihedrals'][i] = adddict['improper']
                break
        else:
            print('added term: ', adddict['improper'])
            copydict['dihedrals'].append(adddict['improper'])

    return copydict

def topol_remove_terms(topoldict,rmvdict):
    copydict = deepcopy(topoldict)
    mapsection = {'pairs':2,'bonds':2,'angles':3,'dihedrals':4}
    for section in ['pairs','bonds','angles','dihedrals']:   
        mapval = mapsection[section]
        for term in rmvdict[section]:
            try:
                searchlist = list(map(itemgetter(*list(range(mapval))),copydict[section]))       # creating list of 
                print(f"removed term {term}")
                del copydict[section][searchlist.index(tuple(term[:mapval]))]
            except ValueError:
                print(f"Couldn't remove term {term} because it is not part of the topology")

    try:
        searchlist = list(map(itemgetter(*list(range(4))),copydict['dihedrals']))
        del copydict['dihedrals'][searchlist.index(tuple(rmvdict['improper'][:4]))]
        print(f"removed term {rmvdict['improper']}")
    except ValueError:
        print(f"Couldn't remove term {rmvdict['improper']} because it is not part of the topology")
    return copydict



        


#%%
## Running the script

input = ConversionRecipe()
input.type = ConversionType.MOVE
input.atom_idx = ['11','9']     #from, to   ## from is the hydrogen that gets moved, to is the heavy atom it forms a bond with
#is_after_homolytic = [[840,842],'m']

toppath = Path("/hits/fast/mbm/hartmaec/workdir/TopologyChanger/example_ala/Ala_delHA_in.top")
ffdir = Path("/hits/fast/mbm/hartmaec/workdir/TopologyChanger/example_ala/amber99sb-star-ildnp.ff")
basedir = Path("/hits/fast/mbm/hartmaec/workdir/TopologyChanger/")
outpath_tmp = Path("/hits/fast/mbm/hartmaec/workdir/TopologyChanger/out/Ala_out_tmp.top")
outpath = Path("/hits/fast/mbm/hartmaec/workdir/TopologyChanger/out/Ala_out.top")
topoldict = read_topol(toppath)
write_topol(topoldict,Path("/hits/fast/mbm/hartmaec/workdir/TopologyChanger/out/Ala_out.top"))

print(topoldict.keys())
print(topoldict['bonds'][1])

if input.type == ConversionType.MOVE:
    print(f"\n----- Starting operation type {input.type} on {input.atom_idx} -----\n")

    heavy_idx = find_heavy(topoldict['bonds'],input.atom_idx[0])
    print(heavy_idx)

    print(f"\n---- Dealing with the 'from' part ----")
    fromGraph = localGraph(topoldict,heavy_idx,ffdir)
    fromGraph.construct_graph()
    print(fromGraph.BondList)
    fromGraph.build_PADs()
    termdict = fromGraph.get_terms_with_atom(input.atom_idx[0])
    termdict['improper'] = []
    fromGraph.remove_terms(termdict)
    print(termdict)
    topoldict = topol_remove_terms(topoldict,termdict)
    print(fromGraph.BondList)
    atom_terms = fromGraph.parameterize_around_atom(heavy_idx)                              
    topoldict = topol_add_terms(topoldict,atom_terms)
    fromGraph.order_lists()
    write_topol(topoldict,outpath_tmp)

    print(f"\n---- Dealing with the 'to' part ----")
    toGraph = localGraph(topoldict,input.atom_idx[1],ffdir)
    toGraph.construct_graph()
    print(toGraph.get_AtomList_property('idx'))
    toGraph.add_Atom(Atom('11'))
    toGraph.add_Bond(['9','11'])
    toGraph.build_PADs()
    print(toGraph.AngleList)
    toGraph.order_lists()
    print(toGraph.get_AtomList_property('idx'))
    print(toGraph.BondList)
    print(toGraph.AngleList)
    toGraph.get_atomprop('atomtype')
    toGraph.get_atomprop('atomname')
    toGraph.get_atomprop('resname')
    print(toGraph.AtomList[7])
    toGraph.AtomList[7].atomtype = 'H1'
    print(toGraph.get_AtomList_property('atomtype'))
    atom_terms = toGraph.parameterize_around_atom(input.atom_idx[1])
    #atom_terms = toGraph.get_terms_with_atom(input.atom_idx[1],center=True,add_function=True)
    print('atom_terms: ',atom_terms)
    atom_terms['improper'] = []  
    print(atom_terms['bonds'])
    print('\n',atom_terms['angles'])
    mapsection = {'bonds':slice(0,1),'angles':slice(1,1),'dihedrals':slice(1,2)} 
    for section in ['bonds','angles','dihedrals']:
        for term in atom_terms[section]:
            if heavy_idx in term[mapsection[section]]:
                atom_terms[section].remove(term)
    print(atom_terms)
    topoldict = topol_add_terms(topoldict,atom_terms)
    toGraph.order_lists()
    atom_terms_H = toGraph.get_terms_with_atom(input.atom_idx[0],add_function=True)
    for section in ['bonds','angles','dihedrals']:
         atom_terms_H[section].clear()
    atom_terms_H['improper'] = []  
    print('atom_terms_H: ',atom_terms_H)
    topoldict = topol_add_terms(topoldict,atom_terms_H)
    toGraph.order_lists()
    termdict = {'bonds':[],'pairs':[],'angles':[],'dihedrals':[]}
    termdict['improper'] = ['7','10','9','14','4','180.0000000','4.6024000','2']
    topoldict = topol_remove_terms(topoldict,termdict)
    write_topol(topoldict,outpath)
    print(toGraph.get_AtomList_property('idx'))
    print(toGraph.get_AtomList_property('atomtype'))
    print(toGraph.AtomList[5].bound_to)
    print([i for i in itertools.combinations(toGraph.AtomList[5].bound_to,2)])
    #fromGraph.add_Atom(Atom('22',atomtype='HC'))
    #fromGraph.add_Bond(['9','22'])






    
    # atom_terms = fromGraph.parameterize_around_radical()

    # print(f"\n---- Changing the Topology ----")
    # print(topoldict['dihedrals'])
    # #topoldict = topol_add_terms(topoldict,atom_terms)
    # topoldict = topol_remove_terms(topoldict,atom_terms)
    # print('\n',topoldict.keys())
    # write_topol(topoldict,outpath)
    # print('\n',topoldict['dihedrals'])

    # print('\n\n--- Wrapping Up ---')
    # print('Atomlist: ',fromGraph.AtomList)
    # print(fromGraph.get_AtomList_property('idx'))
    # print(fromGraph.get_AtomList_property('bound_to'))
    # print(fromGraph.AtomList[0].idx)
    # print('BondList: ',fromGraph.BondList)
    # print('PairList: ',fromGraph.PairList)
    # print('AngleList: ',fromGraph.AngleList)
    # print('DihedralList: ',fromGraph.DihedralList)
    # print('\n')
    # print(fromGraph.get_terms_with_atom(str(22)))
    # print(fromGraph.find_missing_terms())
    # print(atom_terms)



# %%

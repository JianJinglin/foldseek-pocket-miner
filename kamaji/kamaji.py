import sys

from rdkit import Chem

from openbabel import openbabel as ob
from openbabel import pybel
from operator import itemgetter
import numpy
from copy import deepcopy
from biochem import bio, chem
from profile import StructureProfile
from metalclassifier import MetalClassifier
from pdbparser import MultiStatePDB


# TODO recognize L-peptide linking
# TODO add data for BIOLOGICAL MODIFIERS (myristoilation, PLM-palmoilation, ACE, and non-standard, COA, DGA)
#      FAR: farnesyl
# TODO modifiers: BTI (biotinilation)

# TODO phospholypids
#

# TODO add alternate conformation classification (for flexibility analysis?)


# NONSTD AMINO TO RECOGNIZE (double-connection):
# KCX, AHB, OCS, HIC

# PROBLEMATIC STRUCTURES:
# > 5iez

types = ['protein', 'dna', 'rna', 'cofactor',
         'glyco', 'ligand', 'water',
         'metal', 'metal_catalytic', 'salt',
         'modifier', 'additive',
         ]

class Kamaji:
    """ Analyze an OBMol and  a PDB and returns splitted items ready to be
        processed as desired


        TODO: NMR EXMPLE 2l6x.pdb
             - covalent residue!

        REPAIR MECHANISM
              ask the user to re-attach chains that
              are found to be detached or
              residues that are considered ligands...

        ======= WORKING ================

        GUESSING GROUPS
            glycosilation
            lipid acid
            lipid (generic)

        KNOWN CATALYTIC METALS
            Mg, Mn, Fe, Co, Ni, Cu, Zn

        KNOWN SALTS
            Na, Cl, Ca, K, Li, Cl

        KNOWN ADDITIVES
            GOL, SDS (page), LDA

        KNOWN COFACTORS
            Heme, ADP/ATP, GDP/GTP, NADP/NADHI,
            Pyridoxal phosphate, FAD/FADm,

        KNOWN X-RAY ADDITIVES
            Glycerol, SDS, tris-buffer, lauryl dimethylamine-N-Oxide, DMSO,
            sulphate, phosphate, acetate, HEPES, ammonium, ethilen-glycol,
            PEG, ethanol, B-mercapthoethanole, formic acid, pyruvic acid,
            morpholinoethilensulfonate, bicarbonate

        FILTERING
            standard residues
            defined modified residues
            . . . . . . . . . .
            waters
            ions
            cofactors
            additive
            undefined modifiers
            ligands
    """

    def __init__(self, max_short_chain_len = 10, debug=False):
        """
        # XXX ADD: HEADER, COMPND, SOURCE
        # XXX manage peptoids 2SEM
        """
        # any chain with this many residues will be considered a potential ligand
        self._MINRESCOUNT = max_short_chain_len # decapeptide is a ligand
        self._MAXRINGSIZE = 7 # anything >= will be considered macrocycle?
        # self._init_known_types()
        self._init_tools()
        self.graph = {}
        # fingerprints pool
        self._fp_pool = ['fp2', 'fp3', 'fp4', 'maccs']


    def _init_tools(self):
        """ initialize tools and variables to be used later"""
        self._init_known_types()
        # OB tools are initialized here
        self.matcher = ob.OBSmartsPattern()
        self.ob_loader = ob.OBConversion()
        self.ob_loader.SetInAndOutFormats('smi','can')

    def _init_known_types(self):
        """ initialize names and info for known aminoacids,
            nucleic acid residues, cofactors, metals, salts
            and other common PDB ligands
        """
        self.structure_types = {}
               #  'protein': {}, 'nucleic': {}, 'unnown': {}...
        self.structure_groups = {} # clusters, glycosilation sites

        # res_id, numAtoms (max: bound residue count + 1, if terminal)
        self._known_aa = bio.aa_atom_count
        self._known_aa_bonds = bio.aa_bond_count
        self._known_dna = bio.dnabasesAtomCount
        self._known_dna_bonds = bio.dnabasesBondCount
        self._known_rna = bio.rnabasesAtomCount
        self._known_rna_bonds = bio.rnabasesBondCount
        self._known_waters = bio.waterNames
        self._known_metals = chem.metals
        self._known_coordinating_metals = chem.coordinatingMetals
        self._known_salts = chem.salts
        self._known_cofactors = bio.cofactors
        # TODO dead-end...
        # NAG is sugar, PLP should be covalent cofactor,right?
        self._UNDEFINEDMODIF = ['NAG', 'PLP' ] # these are attached to parts of the protein!
        # XXX TOFIX        self._knownResidues = self._known_aa.keys() + self._known_dna.keys() + self._known_rna.keys()
        # SMILES of know bioadditives
        self._known_additives = bio.additives
        # SMARTS of known chemical groups
        self._guessing_groups = chem.simpleorganics

    def sprint(self, text, dest='out'):
        """ stdout stderr"""
        if dest == 'out':
            func = sys.stdout
        elif dest == 'err':
            func = sys.stderr
        func.write(text)

    def set_molecule(self, mol, pdb_info={}):
        """ use a pre-existing multistate PDB object as molecule"""
        self.mol = mol
        self.profile = StructureProfile(self.mol)
        # print("PROFILE HAS NOW", self.profile.get_types())
        #if perceive_bo:
        self.mol.PerceiveBondOrders()
        self.pdb_info = pdb_info
        self._unknown = []
        self.parse_structure()


    def parse_structure(self):
        """
        """
        #self.structure_types = {} -> replaced by self.profile
        # filter residues (protein, dna, rna, water, other)
        print("PARSE RESIDUES...", end=' ')
        self.scan_residues()
        print("[DONE]")
        # search for short chains that could likely be peptide ligands
        print("PARSE SHORT...", end=' ')
        self.scan_short_chain()
        print("[DONE]")
        print("TODO: check for disconnected peptides?")
        # filter all other residues
        self.filter_other()

    def scan_short_chain(self):
        """ scan the graph structure for short chains
            (default: < 10 AA) that could be likely peptide/nucleic
            ligands or modified ones
        """
        types = self.profile.get_types()
        if types == None:
            return
        pool_type = {'protein': 'short_peptide', 'dna':'short_dna', 'rna':'short_rna'}
        for class_, name in pool_type.items():
            if not class_ in types:
                print("SCAN SHORT: class not in pool", class_, types )
                continue
            short_chain = []
            for c in self.profile.get_chains():
                lenght = 0
                residues = self.profile.get_chain(c)
                for rId, info in list(residues.items()):
                    if class_ in info['type']:
                        lenght += 1
                        short_chain.append(rId)
                    if lenght > self._MINRESCOUNT:
                        del short_chain
                        short_chain = []
                        break
                # TODO here is where the short_peptide should be tagged
                # for retrieval
                # TODO use change_data_type(chain_id=c) ??
                for rId in short_chain:
                    self.profile.change_data_type(res_id=rId, type_=name, unique=True)

    def get_short_chain_properties(self, short_chain):
        # XXX this has to be called later in the end...
        """ short chain properties (peptides/nucleic) properties """
        indices = []
        for res_id, data in list(short_chain.items()):
            indices.append((res_id, data['num']))
        sorted(indices, key=itemgetter(1))
        seq = "-".join( [str(x[1]) for x in indices] )
        return { 'sequence' : seq }

    def filter_ligands(self):
        """ scan for potential ligands and known classes"""
        if len(self._unknown) == 0:
            return
        remove = []
        for res_id in self._unknown:
            data = None
            # TODO
            # TODO add here to check if it's a peptide
            # TODO
            lipid = self._is_lipid(res_id)
            if not lipid == False:
                data = {'info': {'class':lipid}}
            else:
                sugar = self._is_glyco_group(res_id)
                if not sugar == False:
                    data = {'info': {'class':'carbohydrate-like'}}
                    # XXX TO DO ImPROVE SuGAR CLASS
                else:
                    data = {'info': {'class':'generic'}}
            if not data == None:
                data['info'].update(self.get_ligand_properties(res_id))
                #info['extra'] = extra
                remove.append(res_id)
                self.profile.add_residue(res_id, 'ligand', data)
        for r in remove:
            self._unknown.remove(r)

        # DEBUG:
        if False:
            for kind, data in list(self.structure_types['ligand'].items()):
                for k in list(data.keys()):
                    #res = self.mol.GetResidue(k)
                    #mol = self._res_to_obmol(res)
                    rObj = self.mol.GetResidue(k)
                    btotal = 0
                    # DEBUG
                    for a in ob.OBResidueAtomIter(rObj):
                        bc =  len([x for x in ob.OBAtomBondIter(a) ])
                        btotal += bc
                    print("RESIDUE TOTAL BONDS",k, btotal)


    def get_ligand_properties(self, res_id):
        """ classifies a generic ligand"""
        mol = self._res_to_obmol(res_id)
        mw = mol.GetMolWt()
        size = 'drug-size'
        if mw < 310:
            size = 'fragment'
        elif mw > 600:
            size = 'large'
        info = {'size': size, 'mw': mw }
        macrocycle_info = self._is_macrocycle(mol)
        if not macrocycle_info is None:
            info['macrocycle'] = macrocycle_info
        return info

    def _is_macrocycle(self, mol):
        """ check if the molecule has a macrocycle ring """
        macro = []
        rings = mol.GetSSSR()
        for r in mol.GetSSSR():
            size = len(r._path)
            if size >= self._MAXRINGSIZE:
                macro.append(size)
        if len(macro):
            return macro
        return None

    def debugwrite(self, fname, reslist=[]):
        c = 0
        for r in reslist:
            #res = self.mol.GetResidue(r)
            mol = self._res_to_obmol(r)
            pybel.Molecule(mol).write('pdb', 'DEBUG_%s_%s.pdb' %(fname, c), overwrite=1)
            c+=1

###     # XXX this one should be moved to a repair object
###     def transformResidues(self):
###         """ """
###         # TODO
###         transformations = biochem.bio.residuenormalization
###         # TYS: sulfotyrosine
###         # phosphotyrosyne
###         # XXX Arsenic
###         # MET C-S-C (1.8) C-Se-C (1.9)

    def filter_other(self):
        """ filter unknown residues into:
            metal, other
        """
        # TODO this needs to be carefully tested to
        #      define the optimal order of application
        #      for each filter. For example PLP (1aam) is a cofactor
        #      and a modifier, which will be classified differently depending
        #      on the order
        if len(self._unknown) == 0:
            return
        # check if the residue is in the Golden Hundred
        # (list of all the top represented ligands
        # in the PDB)
        print("Filter waters")
        self.filter_waters()
        # ions (mono-atomic residues)
        print("Filter ions")
        self.filter_ions()
        print("Filter cofactors")
        self.filter_cofactors()
        print("Filter additives")
        self.filter_additives()
        # at this point, easy groups are all catched
        # remaining ones will be clustered
        print("Cluster unknown...")
        self.cluster_unknown()
        # find ligands (this should be called as last)
        print("Filter ligands...")
        self.filter_ligands()
        # anything that is attached to the protein
        # i.e.: glycosilation, non-std residues
        print("Filter modifiers...")
        self.filter_modifiers()

    def filter_waters(self):
        """ filter waters
            NOTE check here if it is deuterated?
        """
        remove = []
        for res_idx in self._unknown:
            if self._is_water(res_idx):
                remove.append(res_idx)
        for r in remove:
            self._unknown.remove(r)

    def _is_water(self, res_id):
        """ identify water (including deuterated) """
        res = self.mol.GetResidue(res_id)
        atom_elements = [ x.GetAtomicNum() for x in ob.OBResidueAtomIter(res) ]
        atom_elements.sort()
        mono_water = (len(atom_elements)==1 and (8 in atom_elements))
        full_water = (len(atom_elements)==3 and (8 in atom_elements) and (1 in atom_elements) )
        if not mono_water and not full_water:
            return False
        #if (len(atom_elements) == 1) and (not atom_elements == [8]):
        #    return False
        #if (len(atom_elements) == 3) and (not atom_elements == [1,1,8]):
        #    return False
        info = self.get_res_info(res)
        data = {}
        if info['name'] == 'DOD':
            # heavy water
            data['nonstd'] = {'type': 'deuterated', 'standard': 'HOH'}
        #print("WATER IS NOW", res_id, atom_elements, res.GetName(), res.GetNum())
        self.profile.add_residue(res_id, 'water', data)
        return True


    def filter_ions(self):
        """ filter mono-atomic ions and classifies them as
                catalytic
                salt
                metal
        """
        # XXX "[ this function will be updated with statistics on teh PDB distro ]"
        if len(self._unknown) == 0:
            return
        identified = []
        for res_id in self._unknown:
            type_ = None
            res = self.mol.GetResidue(res_id)
            info = self.get_res_info(res)
            #print("SCANINNG", info, res.GetIdx(), res.GetName())
            atoms = [ x for x in ob.OBResidueAtomIter(res) ]
            for x in atoms:
                a = self.mol.GetAtom(x.GetIdx())
                #print("ITER", a.GetIdx(), a.GetAtomicNum())
            #print("AXOMO", atoms[0].GetIdx(), atoms[0].GetAtomicNum())
            if len(atoms) == 1:
                atom = atoms[0]
                #print("ATOM SINGLE", atom.GetIdx())
                if self._is_catalytic_metal(atom): # catalytic
                    type_ = 'metal_catalytic'
                    #print("CATALO")
                # XXX THIS WILL BE CHANGED WITH THE AVaiLABILITY
                # of distributions from the PDB
                elif self._is_salt(atom): # salt
                    type_ = 'salt'
                elif self._is_generic_metal(atom): # generic metal
                    #print("METALO")
                    type_ = 'metal'
            if not (type_ == None):
                self.profile.add_residue(res_id, type_) #, info)
                identified.append(res_id)
        for i in identified:
            self._unknown.remove(i)

    def filter_modifiers(self):
        """ filter anything that is attached to the protein
            as a 'modifier' but was not defined in the
            MODRES keys: glycosilation, very-non-std residues

            NOTE the code assumes that glycosylations and generic
            modifiers are separate groups with no bonds in between.
            This could cause problems in case of fancy glyco-modifiers?
        """
        current_types = self.profile.get_types()
        pool = self._unknown[:]
        for t in ('ligand', 'short_peptide'):
            if not t in current_types:
                continue
            pool += self.profile.get_type(t)
        # for res_id in self._unknown:
        # print("POOL RESIDUES FOR MODIFIER", pool)
        if len(pool) == 0:
            return
        # cache all target atom indices
        target_atoms = []
        for kw in ['protein', 'dna', 'rna']:
            found = self.profile.get_type(kw)
            if found == None:
                continue
            for res in found:
                target_atoms += self._get_atoms_in_res(res)
        target_atoms = set(target_atoms)
        identified = []
        glyco = {}
        # check if residues are attached to the target
        for res_id in pool:
            # print("MODIFIERS PROTECCSING", res_id)
            modifier_data = {}
            res = self.mol.GetResidue(res_id)
            modifier_atom_idx = set(self._get_atoms_bound_to_res(res_id))
            common = target_atoms & modifier_atom_idx
            #print("MODIFI", res_id, common)
            if len(common):
                # print("MODIFIER COMMON PROCESSING", res_id)
                # get residue atom info:
                # atom_info = [self._get_atom_info(x) for x in common]
                # find id of residue attached to the modifier
                #attached_res_idx = self._atom_idx_to_res_idx(common)# [0]
                target_res_idx = self._atom_idx_to_res_idx(common)# [0]
                # find modifier atoms attached to residue
                atoms_bound_to_modified_res = []
                for x in target_res_idx:
                    atoms_bound_to_modified_res += self._get_atoms_bound_to_res(x)
                modified_res_atom_name = [self._get_atom_name(x) for x in common]
                tmp_set_res = set(atoms_bound_to_modified_res)-set(common)
                modifier_atom = list(tmp_set_res & modifier_atom_idx)
                modifier_element = self._obelement_from_atom_idx(modifier_atom)
                modifier_data['group'] = [res_id] + target_res_idx
                modifier_data['target'] = target_res_idx
                modifier_data['target_id'] = [ self.profile.get_residue(x)['id'] for x in target_res_idx]
                modifier_data['modifier_element'] = modifier_element
                # TODO add a check that the modification is done with the oxygen atom?
                if self._is_glyco_group(res_id): # glyco modifier
                    glyco[res_id] = modifier_data['group'] #{attached_res_idx}
                    type_='glyco'
                    # TODO data['info'] = self.getGlycoInfo()
                # TODO add data for BIOLOGICAL MODIFIERS (myristoilation, ACE, and non-standard, COA)
                else: # generic modifier # XXX this should go into filterKnownModifiers?
                    # print("NOT GLYCO!", res_id)
                    modifier_data['info'] = self.get_modifier_info(res_id)
                    type_='modifier'
                # add the modifier information in the profiler
                self.profile.add_residue(res_id,type_, modifier_data)
                # notify the modified residue that it's indeed modified
                for att_idx, att in enumerate(target_res_idx):
                    attached_info = self.profile.get_residue(att)
                    type_ = attached_info['type']
                    if not 'modified' in attached_info:
                        attached_info['modified'] = []
                    #attached_info['modifier_elements'] = self._obelement_from_atom_idx(attaching_atom)
                    attached_info['modifier_elements'] = modifier_element
                    attached_info['modified_atom_name'] = modified_res_atom_name[att_idx]
                    #print("attached mod element", attached_info['modifier_elements'])
                    attached_info['modified'].append(res_id)
                    attached_info['modified'] = list(set(attached_info['modified']))
                    self.profile.add_residue(att, type_, attached_info, replace=True)
                identified.append(res_id)
        # remove known modifier residues
        for i in identified:
            if i in self._unknown:
                self._unknown.remove(i)
        # group glyco modifiers
        if len(glyco):
            self.process_glyco_groups(glyco)

    def _obelement_from_atom_idx(self, idx_list):
        """ return atom symbols (i.e., C, Mg, F...) from a list of indices """
        symbols = []
        for i in idx_list:
            an = self.mol.GetAtom(i).GetAtomicNum()
            symbols.append(ob.GetSymbol(an))
        return symbols

    def process_glyco_groups(self, glyco):
        """ process glycosylation groups and complete polysaccharide structures"""
        # cluster glyco-groups
        print("GLYCO", glyco)
        group_candidates = self._complete_glyco_chains(glyco, self._unknown_clusters)
        found = []
        for group in group_candidates:
            if len(set(group) & set(glyco.keys())) > 0:
                print("GROUP IS", group)
                for g in group:
                    if g in self.profile.data:
                        print("GGGG", g)
                        pp(self.profile.data[g])
                        if 'group' in self.profile.data[g]:
                            current = self.profile.data[g]['group']
                            new = list(set(current+group))
                            self.profile.data[g]['group'] = new
                        else:
                            self.profile.data[g]['group'] = group
                    else:
                        self.profile.add_residue(g, 'glyco', {'group': group})
                found += group
        for f in found:
            try:
                self._unknown.remove(f)
            except:
                pass
        return

    def _complete_glyco_chains(self, glyco, clusters):
        """ group all connected glycosilation residues"""
        glyco_mafia = []
        for res_id in glyco:
            pairs = self._walk_graph( clusters,
                start=res_id,  visited=[] )
            glyco_mafia.append(pairs)
        return glyco_mafia

    def get_modifier_info(self, res_id):
        """ extract information about a modifier residue"""
        res_obj = self._res_to_obmol(res_id)
        res_smiles = self.get_smiles(res_obj)
        res_fp = self._get_mol_fp(res_smiles)
        mw = res_fp['mw']
        if mw < 300:
            size = 'fragment'
        elif 300 < mw < 600:
            size = 'ligand'
        elif mw > 600:
            size = 'large ligand'
        info = {'mw':res_fp['mw'], 'size':size, 'smi':res_fp['smi']}
        return info

    def _get_atom_resnames(self, indices):
        """ return resnames for all atom indices provided"""
        res_names = []
        for i in indices:
            res_obj = self.mol.GetAtom(i).GetResidue()
            try:
                icode = res_obj.getInsertionCode().replace('\x00', '')
            except:
                icode = ""
            res_names.append( '%s%s%s' % (
                res_obj.GetName(), res_obj.GetNumString().strip(), icode ) )
        return res_names

    def _atom_idx_to_res_idx(self, atom_idx_list):
        """ return residue OB indices for all atom OB indices in atom_idx_list"""
        res_idx_list = []
        for i in atom_idx_list:
            res_obj = self.mol.GetAtom(i).GetResidue()
            res_idx_list.append(res_obj.GetIdx())
        return list(set(res_idx_list))

    def _get_atom_name(self, atom_idx):
        """ retrieve atom info from the molecule """
        atom = self.mol.GetAtom(atom_idx)
        res = atom.GetResidue()
        atname = res.GetAtomID(atom)
        return atname.strip()

    def cluster_unknown(self):
        """ generate graph of connected residues
            (i.e. poly-saccharides, n-peptides...)
        """
        if len(self._unknown) == 0:
            return
        clusters = {}
        for r1 in self._unknown:
            clusters[r1] = []
            for r2 in self._unknown:
                if r1 == r2: continue
                if self._are_res_connected(r1, r2):
                    clusters[r1].append(r2)
        self._unknown_clusters = clusters

    def _walk_graph(self, graph,  start, visited = [] ):
        """ recursive function to walk graph
            nodes to find connected clusters
        """
        if not start in visited:
            visited += [start]
        for children in graph[start]:
            if not children in visited:
                visited.append(children)
                new = self._walk_graph( graph, children, visited)
                for n in new:
                    if not n in visited:
                        visited.append(n)
        return visited

    def _are_res_connected(self, res1, res2):
        """ res1, res2:
            check if two residues are connected by
            at least a bond between two atoms
        """
        atoms_bound_to_r1 = set( self._get_atoms_bound_to_res(res1) )
        atoms_in_r2 = set( self._get_atoms_in_res(res2) )
        return len( atoms_bound_to_r1 & atoms_in_r2 ) > 0

    def _get_atoms_in_res(self, res_id):
        """ convert a list of residues to the list of indices
            of all atoms contained in the residues
        """
        return [ x.GetIdx() for x in self.get_residue_atoms(res_id) ]

    def _get_atoms_bound_to_res(self, res_idx, ignore_metals=True):
        """ find all atoms with which atoms in the
            residue establish bonds
        """
        # TODO this does not work for metals!
        bound = []
        atom_list = self.get_residue_atoms(res_idx)
        for a in atom_list:
            symbol = ob.GetSymbol(a.GetAtomicNum())
            if symbol in self._known_metals:
                continue
            bound += [ b.GetBeginAtomIdx() for b in ob.OBAtomBondIter(a) ]
            bound += [ b.GetEndAtomIdx() for b in ob.OBAtomBondIter(a) ]
        bound = set(bound)
        if ignore_metals:
            out = []
            for i in bound:
                symbol = ob.GetSymbol(self.mol.GetAtom(i).GetAtomicNum())
                if symbol in self._known_metals:
                    continue
                out.append(i)
            return out
        else:
            return bound

    def _is_protein_res(self, res_obj):
        """ check that a residue belongs to protein and collect info

            NOTE: actual residues that are tagged as HET in the protein entry
            will be considered as modifiers (e.g., A:CYS708 in 2ibn)
        """
        res_id = res_obj.GetIdx()
        rinfo = self.get_res_info(res_obj)
        rchain = rinfo['chain']
        rname = rinfo['name']
        rnum = int(rinfo['num'])

        if not rname in self._known_aa:
            return False
        # check if it's a residue playing another role
        for a in ob.OBResidueAtomIter(res_obj):
            if res_obj.IsHetAtom(a):
                return False
        data = {}
        # mutated residue
        mut_info = self.pdb_info.get('seqAdv',[])
        if not len(mut_info) == 0:
            if rchain in mut_info:
                if rnum in mut_info[rchain]:
                    if not rname == mut_info[rchain][rnum]['res_name']:
                        print("Warning: mismatch between PDB header info and structural data")
                    data['mutation'] = {'org': mut_info[rchain][rnum]['std_res_name'],
                            'description': mut_info[rchain][rnum]['mod_desc'],
                            }
        missing_atom, missing_bond = self._check_missing_atom_bond(res_obj, type_='protein')
        # a tolerance of 1 is added to atoms to account for terminal amino acids
        if (missing_atom > 1) or (missing_bond > 0) :
            data['missing'] = {'fixed': False, 'atoms':missing_atom-1, 'bonds':missing_bond}
        self.profile.add_residue(res_id, 'protein', data)
        return True

    def _is_dna_res(self, res_obj):
        """ check that a residue belongs to DNA and collect info
            NOTE: atom and bond counts take into account the possible
                  absence of a phosphate group (+/4 atoms and +/-3 bonds)
        """
        res_id = res_obj.GetIdx()
        data = {}
        info = self.get_res_info(res_obj)
        name = info['name']
        if not name in self._known_dna:
            return False
        missing_atom, missing_bond = self._check_missing_atom_bond(res_obj, type_='dna')
        # a tolerance of 4 atoms and 3 bonds is added for terminal phosphates
        if (missing_atom > 4) or (missing_bond > 3):
            data.update({ 'missing': {'fixed':False, 'atoms':missing_atom-4, 'bonds':missing_bond}})
        self.profile.add_residue(res_id, 'dna', data)
        return True

    def _is_rna_res(self, res_obj):
        """ check that a residue belongs to DNA and collect info
            NOTE: atom and bond counts take into account the possible
                  absence of a phosphate group (+/4 atoms and +/-3 bonds)
        """
        res_id = res_obj.GetIdx()
        data = {}
        info = self.get_res_info(res_obj)
        name = info['name']
        if not name in self._known_rna:
            return False
        missing_atom, missing_bond = self._check_missing_atom_bond(res_obj, type_='rna')
        # a tolerance of 4 atoms and 3 bonds is added for terminal phosphates
        # TODO fix the potentially missing phosphate (1yy0)
        if (missing_atom > 4) or (missing_bond > 3):
            data.update({ 'missing': {'fixed':False, 'atoms':missing_atom-4, 'bonds':missing_bond}})
        self.profile.add_residue(res_id, 'rna', data)
        return True

    def _check_missing_atom_bond(self, res_obj, type_):
        """ checks for missing atoms and bonds in standard residues and nucleobases """
        data_pool = { 'protein': (self._known_aa, self._known_aa_bonds),
                    'dna': (self._known_dna, self._known_dna_bonds),
                    'rna': (self._known_rna, self._known_rna_bonds)}
        # gather info about the current residue
        atom_count, bond_count = data_pool[type_]
        rinfo = self.get_res_info(res_obj)
        rchain = rinfo['chain']
        rname = rinfo['name']
        rnum = int(rinfo['num'])
        atom_current = self.count_atom_in_res(res_obj)
        bond_current = self.count_bonds_in_res(res_obj)
        # count missing parts
        missing_atom = atom_count[rname] - atom_current
        missing_bond = bond_count[rname] - bond_current
        return missing_atom, missing_bond

    def _is_non_std_res(self, res_obj):
        """ check if a residue is in the set of
            MODRES residues defined in the header
            of the PDB
        """
        # TODO this may be fragile, check with ProDy info
        res_id = res_obj.GetIdx()
        if not 'modRes' in self.pdb_info:
            return
        found_mod = self.pdb_info['modRes']
        info=self.get_res_info(res_obj)
        chain = info['chain']
        name = info['name']
        num = info['num']
        icode = info['icode']
        data = {}
        if not chain in list(found_mod.keys()):
            return False
        for mod_res in found_mod[chain]:
            if not ( mod_res['res_name'] == name ):
                continue
            if not ( mod_res['res_num'] == num ):
                continue
            if not ( mod_res['iCode'] == icode ):
                continue
            data = {'nonstandard' : { 'type' : mod_res['mod_description'],
                                 'standard' : mod_res['std_res_name'] }}
            break
        if data == {}:
            return False
        std = data['nonstandard']['standard']
        if std in self._known_aa:
            type_ = 'protein'
        elif std in self._known_dna:
            type_ = 'dna'
        elif std in self._known_rna:
            type_ = 'rna'
        else:
            type_ = 'unknown'
        self.profile.add_residue(res_id, type_, data)
        return True

    def get_residue_atoms(self, res_idx):
        """ return all atoms in a residue"""
        res = self.mol.GetResidue(res_idx)
        return [ x for x in ob.OBResidueAtomIter(res) ]

    def _is_glyco_group(self, res_id):
        """ check if a residue is
            a glycosilating group
        """
        pattern = self._guessing_groups['glycosilation']
        resMol = self._res_to_obmol(res_id)
        result = self.SMARTSMatcher(resMol, pattern)
        #print("PATTERN", pattern, result)
        # XXX TODO ADD A MW ratio?
        if result:
            # TODO classify better the sugar
            return True #'carbohydrate'
        return False

    def _is_lipid(self, res_idx):
        """ check if a residue is a lipid [EXPAND DEFINITION A BIT!]"""
        res = self._res_to_obmol(res_idx)
        patterns = self._guessing_groups['lipids']
        smi = self.get_smiles(res)
        for pair in patterns:
            name, p = list(pair.items())[0]
            result = self.SMARTSMatcher(res, p)
            if result:
                return name
        return False

    def SMARTSMatcher(self, mol, pattern):
        # XXX convert this to OBSmarts
        """ return the matching indices for
            pattern in mol.
            Mol can be a OBMol or a OBResidue
        """
        #matcher = ob.OBSmartsPattern()
        self.matcher.Init(pattern)
        # XXX TODO change this to OB
        #return matcher.Match(mol)
        matcher = pybel.Smarts(pattern)
        return matcher.findall(pybel.Molecule(mol))


    def _find_most_similar(self, res_id, known, score_cutoff=0.5,
            use_chem_formula=False):
        """ identify the closest guess for res_id in a given known class
            using a multi-fingerprint tanimoto similarity comparison
        """
        #from pprint import pprint
        res_obj = self._res_to_obmol(res_id)
        res_smiles = self.get_smiles(res_obj)
        res_fp = self._get_mol_fp(res_smiles)
        res_chem_formula = self._get_chem_formula(mol=res_obj)
        best_score = -1
        best_chem_formula = 1000
        best_name = None
        for guess_id, data in list(known.items()):
            #print("GUIESS", guess_id)
            name = data['name']
            smi = data['smi']
            guess_fp = self._get_mol_fp(smi)
            result = self._calc_similarity(res_fp, guess_fp)
            score = result['score']
            guess_chem_formula = self._get_chem_formula(smi=smi)
            score_chem_formula = self._calc_chem_formula_difference(guess_chem_formula, res_chem_formula)
            #print("SCORES:", score_chem_formula, "RES:", result)
            if score > best_score:
                best_score = score
                best_id = name
                best_match = result
                best_smiles = smi
                best_info = data
            best_chem_formula = min(score_chem_formula, best_chem_formula)
        if best_score < score_cutoff:
            return None
        if use_chem_formula:
            if (best_chem_formula > 2) and best_score < 0.85:
                #print("\n\n\n\nFAILURE!!!!\n\n\n")
                return None
        accuracy = 'low'
        if best_score >=0.7:
            accuracy = 'high'
        guess_info = {'estimate_accuracy':accuracy, 'estimate_score':best_score,
                'estimate':best_id, 'match':best_match, 'smiles':best_smiles,
                'name':best_name, 'known':best_info, 'smiles_query': res_smiles}
        #print("\n\nWE WILL RETURN GUESS INFO", guess_info)
        return guess_info

    def _get_chem_formula(self, mol=None, smi=None):
        """ perform elemental analysis to check if  """
        if not smi == None:
            mol = ob.OBMol()
            self.ob_loader.ReadString(mol, smi)
        elements = {}
        for a in ob.OBMolAtomIter(mol):
            enum = a.GetAtomicNum()
            elements.setdefault(enum, 0)
            elements[enum] += 1
        return elements

    def _calc_chem_formula_difference(self, elements1, elements2):
        """ score based on brute formula analysis, the lower the better
            one extra element type  : +1
            one extra element count : +1
        """
        score = 0
        for e1, count1 in list(elements1.items()):
            if not e1 in elements2:
                score += 1
            else:
                score += abs(elements1[e1] - elements2[e1])
        for e2, count2 in list(elements2.items()):
            if not e2 in elements1:
                score += elements2[e2]
        return score


    def filter_additives(self):
        """ filter for known solvent/additives

            this should rely on the top 100 most represented
            ligands in the PDB
        """
        # this function may be streamlined
        if len(self._unknown) == 0:
            return
        #cutoffFP = 0.8
        #cutoffMWratio = 0.5
        score_cutoff = 0.65
        found = []
        for res_id in self._unknown:
            closest = self._find_most_similar(res_id, self._known_additives, score_cutoff)
            if closest == None:
                continue
            res = self.mol.GetResidue(res_id)
            data = {'info': closest}
            self.profile.add_residue(res_id, 'additive', data)
            found.append(res_id)
        for f in found:
            self._unknown.remove(f)

    def get_smiles(self, mol):
        """ generate canonical smiles of OBMol mol"""
        return self.ob_loader.WriteString(mol).strip()

    def _get_mol_fp(self, smi):
        """ """
        #mol = ob.OBMol()
        #self.ob_loader.ReadString(mol, smi)
        mol = pybel.readstring('can', smi)
        data = {'smi' : smi, 'mw': mol.molwt, 'fp' : {} }
        for fp in self._fp_pool:
            data['fp'][fp] = mol.calcfp(fp)
        return data

#    # OB VERSION... utterly complicated
#    def _calc_fp(self, mol, fp_type):
#        """ calculate molecular descriptors """
#        fp = ob.VectorUnsignedInt()
#        ...forgetaboutit...



    def _calc_similarity(self, mol1_data, mol2_data):
        """ calculate score between two sets of fingerprint"""
        mw1 = mol1_data['mw']
        mw2 = mol2_data['mw']
        if (mw1 == 0) or (mw2 == 0):
            return {'score': 0.0, 'score_avg': 0.0,
                'tanimoto_avg' : 0.0, 'mwratio' : 0.0,
                'tanimoto': 0.0, 'tanimoto_best':0.0}
        similarity = {}
        avg = 0
        tanimoto_best = -1
        for fp_type in list(mol1_data['fp'].keys()):
            fp1 = mol1_data['fp'][fp_type]
            fp2 = mol2_data['fp'][fp_type]
            tanimoto = fp1 | fp2
            similarity[fp_type] = tanimoto
            avg += tanimoto
            tanimoto_best = max(tanimoto_best, tanimoto)
        avg = avg / len(mol1_data['fp'])
        try:
            mw_ratio = 1 - ( abs(mw1 - mw2) / mw1)
        except ZeroDivisionError:
            print(mw1, mw2)
            mw_ratio = 0.01
        score_avg = avg * mw_ratio
        #print "BEST mw_ratio %2.2f %2.2f" % (tanimoto_best, mw_ratio)
        score = tanimoto_best * mw_ratio
        return {'score': score, 'score_avg': score_avg,
            'tanimoto_avg' : avg, 'mw_ratio' : mw_ratio,
            'tanimoto': similarity, 'tanimoto_best':tanimoto_best}

    def _res_to_obmol(self, res_id, forcesingle=False):
        """ create a separate molecule from a
            residue with res_id within a molecule
        """
        residue = self.mol.GetResidue(res_id)
        mol = ob.OBMol()
        mol.BeginModify()
        table = {}
        bonds = []
        for a in ob.OBResidueAtomIter(residue):
            new = mol.NewAtom()
            new.SetAtomicNum( a.GetAtomicNum() )
            new.SetVector( a.GetVector() )
            table[a.GetIdx()] = new.GetIdx()
        miss = 0
        for a in ob.OBResidueAtomIter(residue):
            for b in [ x for x in ob.OBAtomBondIter(a) ]:
                begin = b.GetBeginAtomIdx()
                end = b.GetEndAtomIdx()
                order = b.GetBondOrder()
                if forcesingle:
                    order = 1
                try:
                    mol.AddBond( table[begin], table[end], order)
                except:
                    miss+=1
                    pass
        #print("*** MISSING BONDS *** ", miss)
        mol.EndModify()
        mol.ConnectTheDots()
        mol.PerceiveBondOrders()
        for a in ob.OBMolAtomIter(mol):
            ob.OBAtomAssignTypicalImplicitHydrogens(a)
        #self.ob_loader.SetOutFormat("can")
        #print(self.ob_loader.WriteString(mol))
        return mol

    def blendmol(self, mol):
        """ """
        mol.write('pdb', filename='x.pdb', overwrite=True)
        mol = next(pybel.readfile('pdb', 'x.pdb'))
        mol.OBMol.ConnectTheDots()
        return mol

    def _is_in_chem_class(self, atom, chem_class):
        """ generic function to check if the element of
            an atom belongs to a chemical class,
            (list of element symbols)
        """
        eNumber = atom.GetAtomicNum()
        eSymbol = ob.GetSymbol(eNumber)
        #print("CLESS", eNumber, eSymbol,chem_class)
        return eSymbol in chem_class

    def _is_catalytic_metal(self, atom):
        """ check if an atom is in the list of metals
            known to have a catalyric role
        """
        klass =  list(self._known_coordinating_metals.keys())
        if self._is_in_chem_class(atom, klass):
            # XXX TODO check the coordination geometry
            return True
        return False

    def _is_generic_metal(self, atom):
        """ check (bool) if an atom is a metal"""
        return self._is_in_chem_class(atom, self._known_metals)

    def _is_salt(self, atom):
        """ check if an atom is a salt ion"""
        return self._is_in_chem_class(atom, self._known_salts)

    def _is_cofactor(self, res_id, res_name, smarts=True):
        """ check if a residue is a known cofactor
            (name, for now, SMARTS later?)
        """
        score_cutoff = 0.5 # minimum score to be recognized as cofactor
        known = self._known_cofactors.get(res_name, None)
        if not known == None and not smarts:
            return (res_name, known)
        #print "RESNAME", res_name
        closest = self._find_most_similar(res_id, self._known_cofactors, score_cutoff, use_chem_formula=True)
        if closest == None:
            return None
        return closest

        """
        for cofactorId, data in self._known_cofactors.items():
            name = data['name']
            smi = data['smi']
            cofactorFP = self._get_mol_fp(smi)
            result = self._calc_similarity(resFP, cofactorFP)
            score = result['score_avg']
            if score > bestScore:
                bestScore = score
                bestId = cofactorId
                bestMatch = result
                bestSMI = smi
                bestName = name
                bestData = data
        if bestScore < score_cutoff:
            return None
        accuracy = 'low'
        if bestScore >=0.7:
            accuracy = 'high'
        extra = {'estimate_accuracy':accuracy, 'estimate_score':bestScore,
                'estimate':bestId, 'match':bestMatch, 'smi':bestSMI,
                'name':bestName, 'known':bestData}
        return extra
        """

# OBSOLETE
#    def _isADSupported(self, atom):
#        """ check if is is an AutoDock supported element"""
#        return self._is_in_chem_class(atom, self._autodockSupportedElements)


    def filter_cofactors(self):
        """ filter and identifies known co-factors """
        if len(self._unknown) == 0:
            return
        identified = []
        for res_id in self._unknown:
            res = self.mol.GetResidue(res_id)
            info = self.get_res_info(res)
            found = self._is_cofactor(res_id, info['name'])
            if not found == None:
                data = {'info': found}
                identified.append(res_id)
                self.profile.add_residue(res_id, 'cofactor', data)
        for i in identified:
            self._unknown.remove(i)


    def scan_residues(self):
        """ subdivide residues by type and populate
            the structureType dictionary:
                protein
                nucleic
                other (anything else)

            for each recognized residue from biopolymers
            the number of atoms found is checked with expected
            values plus a tolerance number:
                - aminoacids : +/-1 for terminal residue
                - nucleic  : +/-4 for 3'/5' terminal (PO3-O)
        """
        for res_obj in ob.OBResidueIter(self.mol):
            #print("SCANNING", res_obj.GetName(),res_obj.GetNum())
            if self._is_protein_res(res_obj):
                continue
            elif self._is_dna_res(res_obj):
                continue
            elif self._is_rna_res(res_obj):
                continue
            elif self._is_non_std_res(res_obj):
                continue
            else:
                self._unknown.append( res_obj.GetIdx() )
        return

    def get_res_info(self, res):
        """ retrieve info about a residue"""
        #if not res_idx == None:
        #    res = self.mol.GetResidue(res_idx)
        name = res.GetName().strip()
        num = int(res.GetNumString())
        chain = res.GetChain().strip()
        try:
            icode = res.GetInsertionCode().replace('\x00', '')
        except:
            #print "Warning: no insertion code"
            icode =''
        #return chain, name, num
        return {'chain':chain, 'name':name, 'num':num, 'icode':icode}

    def count_atom_in_res(self, res=None, heavy_only=True):
        """ count number of atoms in a residue
        """
        #if not res_idx == None:
        #    res = self.mol.GetResidue(res_idx)
        c = 0
        for a in ob.OBResidueAtomIter(res):
            if (a.GetAtomicNum() == 1) and heavy_only:
                continue
            c+= 1
        return c

    def count_bonds_in_res(self, res=None, heavy_only=True):
        """ count number of bonds in a residue
        """
        #if not res_idx == None:
        #    res = self.mol.GetResidue(res_idx)
        c = 0
        at_idx_list = []
        bond_idx_list = []
        bond_list = []
        for a in ob.OBResidueAtomIter(res):
            if (a.GetAtomicNum() == 1) and heavy_only:
                continue
            at_idx_list.append(a.GetIdx())
            for b in ob.OBAtomBondIter(a):
                if (b.GetBeginAtom().GetAtomicNum() == 1) and heavy_only:
                    continue
                if (b.GetEndAtom().GetAtomicNum() == 1) and heavy_only:
                    continue
                bond_list.append(b)
        for b in bond_list:
            if b.GetBeginAtomIdx() in at_idx_list:
                if b.GetEndAtomIdx() in at_idx_list:
                    bond_idx_list.append(b.GetIdx())
        return len(set(bond_idx_list))



    def get_residue(self, name=None, num=None, chain=None):
        """ retrieve the residue(s) matching the
            name-num-chain combination)
        """
        found = []
        for res in ob.OBResidueIter(self.mol):
            rname = res.GetName()
            rnum = int(res.GetNumString())
            rchain = res.GetChain()
            if (not name == None) and (not rname == name):
                continue
            if (not num == None) and (not rnum == num):
                continue
            if (not chain == None) and (not rchain == chain):
                continue
            #if (rname == name) or (name == None):
            #    if (rnum == num) or (num == None):
            #        if (rchain == chain) or (chain == None):
            #           found.append(res.GetIdx())
            found.append(res.GetIdx())
        return found


    def fixWater(self):
        """ fix water molecules:
                - add hydrogens to lone oxygens
                - convert deuterium to hydrogen
        """
        # XXX to be written
        pass



if __name__ == '__main__':
    import sys
    import os
    from pprint import pprint as pp

    def res_info(res_idx):
        rObj = K.mol.GetResidue(res_idx)
        info = K.get_res_info(rObj)
        profile = c.get_residue(res_idx)
        print(info, profile)

    def getLines(filename, doStrip = False, removeEmpty=False):
        """ """
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()
        if doStrip:
            lines = list(map(strip,lines))
        if removeEmpty:
            #lines = removeEmptyLines(lines)
            lines = [ l for l in lines if l.strip() ]
        return lines



    def printExpand(data, indent=0, exclude=[]):
        """ print nested data

         [A]
          |
          --[ modifier]
          |      |
          |    [ 57 ]
          |      |
          |      |
        """
        if indent > 0:
            #spacer = "" + ("%s|" % (" " *indent))  * (indent+1)
            spacer = "" + ( "%s|" % ("    " * indent) )
            curve = "" + ( "%s'" % ("    " * indent) )
        else:
            spacer = ""
            curve = ""
        print("\n"+spacer,"\n"+spacer, end=' ')
        #print "\n"+curve,

        #print "SPACER[%s]" % spacer
        if isinstance(data, dict):
            for k in sorted(data.keys()):
                if k in exclude:
                    continue
                v = data[k]
                print("\n%s--[ %s ]" % (curve, k), end=' ')
                if isinstance(v, dict) or isinstance(v,list):
                    printExpand(v, indent+1, exclude)
                else:
                    print("--> ( %s )" % v, end=' ')
        elif isinstance(data, list):
            for v in data:
                if v in exclude:
                    continue
                if isinstance(v, dict) or isinstance(v,list):
                    printExpand(v, indent+1, exclude)
                else:
                    #print "(",v,")"
                    print("%s" % v, end=' ')
                    #print "%s" % (" "*indent), v,


    missing_res = 0
    # load molecule and perform multiplicity analysis
    molfile = sys.argv[1]
    name, ext = os.path.splitext(molfile)
    pdb_id = os.path.basename(name)
    ext = ext[1:].lower()
    if ext == "pdb":
        rawpdb = getLines(molfile)
        #print "[ %d lines]" % len(rawpdb)
        mmol = MultiStatePDB(rawpdb)
        #print "Molecule loaded has %d model(s)" % mmol.modelCount
        if mmol.model_count == 0:
            #print "[ no models found in the structure ]"
            exit(0)
        #print "Molecule loaded has %d residue(s) with alternate locations" % len(mmol.getAltResidues)
        print("Chains found:", mmol.chains)
        mmol.set_state(model=0, alt_mode='a') # model = 0, altMode='A', altResidueMode={ 'B:THR276' : 'A'}
        mol = mmol.get_structure()
        info = mmol.pdb_info
        for chain, miss in list(mmol.pdb_info['missingRes'].items()):
            missing_res += len(miss)
    else:
        parser = ob.OBConversion()
        parser.SetInFormat(ext)
        mol = ob.OBMol()
        parser.ReadFile(mol, molfile)
        info={}
    K = Kamaji()
    # initialize Kamaji and perform structure analysis
    K.set_molecule(mol, pdb_info=info)
    c = K.profile


    extraction = True

    if extraction == True:
        sel_list = []
        import prody
        pmol = prody.parsePDB(sys.argv[1])
        sel_list = []
        for res_id, res_info in c.get_type('protein').items():
            sel_list.append('(resnum {num} and chid {chain})'.format(**res_info))
        if 'cofactor' in c.get_types():
            for co_id, co_info in c.get_type('cofactor').items():
                sel_list.append('(resnum {num} chid {chain})'.format(**co_info))
        # print("SEL STRING", sel_list)
        sel_string = " or ".join(sel_list)
        out_sel = pmol.select(sel_string)
        print("SAVE: %s_target.pdb"% name)
        prody.writePDB('%s_target.pdb' % name, out_sel)

        mod_res_list =[]
        print("TYPES IS", c.get_types())
        print("CCC", c.get_type("modifier"))
        if 'modifier' in c.get_types():
            for mod_id, mod_data in c.get_type('modifier').items():
                target = mod_data['target']
                for t in target:
                    r = c.get_residue(t)
                    mod_id = r['modified'][0]
                    mod = c.get_residue(mod_id)
                    mod_res_list.append("MOD>> %s,%s,%s:%s%d,%s" % (r['name'],pdb_id, r['chain'], r['name'], r['num'], mod['name']))
                    print(mod_res_list[-1])

        sys.exit(0)

    print("Chains found          : %d" % len(c.get_chains()))
    print("Residues found        : %d" % len(list(c.data.keys())))
    print("Known missing residues: %d" % missing_res)
    details = ['short_peptide']
    types = c.get_types()
    print("\n\nFOUND CLASSES")
    for t in types:
        print("- %s" % t)
###     if ('metal_catalytic' in types) or ('cofactor' in types):
###         metal = MetalClassifier(c)
###         print "========================="
###         print "Parsing metals"
###         if 'metal_catalytic' in types:
###             ions = c.get_type('metal_catalytic')
###             print "  ions:  %d" % (len(ions.keys()))
###             metal.scanMetalIon()
###         if 'cofactor' in types:
###             metal_cat = c.get_type('cofactor')
###             print "  potential cofactor-ions:  %d" % (len(metal_cat.keys()))
###             metal.scanMetalCofactor()
###         metal.processAll()
###         print "=========================\n"

    skippo = ["protein", "water"]

    for t in c.get_types():
        if t in skippo:
            continue
        print("[[[[[[[[[[ info %s ]]]]]]]]]]" % t)
        pp(c.get_type(t))

    print("=========================")
    print("FOUND TAGS")
    print(c.get_tags())
    print("=========================")
    tags = c.get_tags()
    if tags == None:
        sys.exit()
    if 'missing' in c.get_tags():
        print("MISSING:")
        mi = c.get_tag('missing')
        for k,v in list(mi.items()):
            #pp(v)
            print('\t residue %s [%d]' % (v['id'],k), 'missing atoms:', v['missing']['atoms'], 'bonds', v['missing']['bonds'])
    if 'modified' in c.get_tags():
        #print("MODIFIED:")
        mo = c.get_tag('modified')
        for k,v in list(mo.items()):
            print('\t residue',k,v['id'], "is modified by", end=' ')
            modifier = v['modified']
            for momo in modifier:
                momo_info = c.get_residue(momo)
                group = [str(x) for x in momo_info['group'] if not x == k]
                print("%s" % ( "-".join(group)))




    sys.exit()

    #s = k.structure_types
    #print "NOTES: mutate Selenium, Arsenic-Cys..."
    #print "NOTES: check if a modifier is some sort of AA to tag it for normalizing mutation to a natural AA"

    print("ADD ARSENIC 4fsd")
    print("To process")
    print("- missing atoms")
    print("- recognize coordinating metals")
    print("- cofactors")
    print("- non-std residues")
    print("- ADD biotinilation")
    sys.exit()


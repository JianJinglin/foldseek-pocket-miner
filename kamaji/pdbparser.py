from openbabel import openbabel as ob



class MultiStatePDB:
    """ 
        deals with MODEL and alternate residue locations...
    
    """
    def __init__(self, text):
        """ """
        self.text = text
        self.init_vars()
        self.parse_lines()
        self.parse_header()
        self.parse_atoms() 

    def parse_header(self):
        """ 
          PDB FILE DEFINED PROPERTIES
            Herarchical structure tree
            Modified residues
            Experimental data
            Missing atoms
            Symmetry
            Biounit
            Hetero list
        """
        self.init_pdb_header_parser()
        self.parse_pdb_info()

    def init_pdb_header_parser(self):
        """ initialize mapping betwen PDB info and parsers"""

        self._type_to_func = { 
                'missingRes' : self._parse_missing_residues,
                'missingAtoms' : self._parse_missing_atoms,
                'symmetry' : self._parse_symmetry,
                'biounit' : self._parse_biounit,
                'hetList': self._parse_het_list,
                'expData' : self._parse_exp_data,
                'modRes' : self._parse_mod_res,
                'ssBond' : self._parse_disulf_bond,
                'seqAdv': self._parse_seq_adv,
              }

        self.remark_to_type = { 
                'REMARK 465' : 'missingRes',
                'REMARK 470' : 'missingAtoms',
                'REMARK 290' : 'symmetry',
                'REMARK 350' : 'biounit',
                'HET   '    : 'hetList',
                'EXPDTA'     : 'expData',
                'MODRES' : 'modRes',
                'SSBOND': 'ssBond',
                'SEQADV': 'seqAdv',
                }

        self.type_to_remark = { 
                'missingRes' : 'REMARK 465',
                'missingAtoms' :'REMARK 470',
                'symmetry' : 'REMARK 290',
                'biounit' : 'REMARK 350',
                'hetList' : 'HET   ',
                'expData': 'EXPDTA',
                'modRes': 'MODRES',
                'ssBond': 'SSBOND',
                'seqAdv' :'SEQADV',
                  }

    def parse_pdb_info(self):
        """ parse all remarks in the PDB"""
        self.pdb_info = {}
        for kw, v in list(self.remark_to_type.items()):
            self.pdb_info[v] = []
        for l in self.get_header():
            value = None
            if l[0:6] == 'REMARK':
                value = l[0:10]
        for l in self.get_header():
            value = None
            if l[0:6] == 'REMARK':
                value = l[0:10]
                if not value in list(self.remark_to_type.keys()): 
                    value = None
                    continue
            elif l[0:6] == 'HET   ':
                value = l[0:6]
            elif l[0:6] =='EXPDTA':
                value = l[0:6]
            elif l[0:6] == 'MODRES':
                value = l[0:6]
            elif l[0:6] == 'SEQADV':
                value = l[0:6]
            #elif l.startswith('ATOM') or l.startswith('HETATM'):
            #    pass
            if not value == None:
                #print "VALUE", value
                kw = self.remark_to_type[value]
                buff = self.pdb_info[kw].append(l)
        
        for kw, data in list(self.pdb_info.items()):   
            self.pdb_info[kw] = self._type_to_func[kw](data)

    def _parse_missing_residues(self, data):
        """ parse the missing residues entry"""
        new = {}
        kw = self.type_to_remark['missingRes']
        pattern = 'REMARK 465   M RES C SSSEQI'
        inside = False
        for l in data:
            if inside:
                #raw = l.split(kw, 1)[1]
                raw = l[15:]
                #print "RAW", l
                res, chain, seq_id = raw.split() # this must be always len() == 3?!?
                if not chain in list(new.keys()):
                    new[chain] = []
                new[chain].append( (res, seq_id) )
            if pattern in l:
                inside = True
        return new

    def _parse_exp_data(self, data):
        """
            return the experimental method
            NOTE: some intelligence would be useful here...
        """
        kw = 'EXPDTA'
        method = 'other'
        known_methods = { 'nmr': 'nmr',
                    'electron crystallography' : 'ec',
                    'x-ray diffraction': 'xray',
                    }
        raw = []

        for l in data:
            raw.append( l.split(kw, 1)[1].strip() )
        raw_string = " | ".join(raw)
        for k,v in list(known_methods.items()):
            if k in raw_string.lower():
                method = v
                break
        return  { 'raw': raw, 'method': v }

    def _parse_missing_atoms(self, data):
        """ """ 
        return data

    def _parse_disulf_bond(self, data):
        """ ,
        [ SOURCE: http://www.wwpdb.org/documentation/format33/sect6.html#SSBOND ]
        COLUMNS        DATA  TYPE     FIELD            DEFINITION
        --------------------------------------------------------------------------------
         1 -  6        Record name    "SSBOND"
         8 - 10        Integer        serNum           Serial number.
        12 - 14        LString(3)     "CYS"            Residue name.
        16             Character      chainID1         Chain identifier.
        18 - 21        Integer        seqNum1          Residue sequence number.
        22             AChar          icode1           Insertion code.
        26 - 28        LString(3)     "CYS"            Residue name.
        30             Character      chainID2         Chain identifier.
        32 - 35        Integer        seqNum2          Residue sequence number.
        36             AChar          icode2           Insertion code.
        60 - 65        SymOP          sym1             Symmetry operator for residue 1.
        67 - 72        SymOP          sym2             Symmetry operator for residue 2.
        74 - 78        Real(5.2)      Length           Disulfide bond distance
        """
        new = {}
        return new
        # see _parse_mod_res for adaptation




    def _parse_mod_res(self, data):
        """ processes modified residues
        [ SOURCE: http://www.wwpdb.org/documentation/format23/sect3.html#MODRES ]

        COLUMNS    DATA TYPE        FIELD         DEFINITION
        ----------------------------------------------------
         1 - 6     Record name      "MODRES"
         8 - 11    IDcode           idCode     ID code of this entry.
        13 - 15    Residue name     res_name    Residue name used in this entry.
        17         Character        chainID    Chain identifier.
        19 - 22    Integer          seqNum     Sequence number.
        23         AChar            iCode      Insertion code.
        25 - 27    Residue name     stdRes     Standard residue name.
        30 - 70    String           comment    Description of the residue
        """
        kw = self.type_to_remark['modRes']
        new = {}
        for l in data:
            pdb_id = l[7:11].strip()
            res_name = l[12:15].strip()
            chain_id = l[16] #.strip()
            res_num = int(l[18:22].strip())
            iCode = l[22].strip()
            std_res_name = l[24:27].strip()
            mod_description = l[29:70].strip()
            chain = new.setdefault(chain_id, [])
            chain.append( { 'res_name': res_name,
                                   'chain_id': chain_id,
                                   'res_num' : res_num,
                                   'iCode'  : iCode,
                                   'std_res_name': std_res_name,
                                   'mod_description': mod_description,
                                 })
        return new

    def _parse_seq_adv(self, data, ignore_deletion=True):
        """ The SEQADV record identifies conflicts between sequence 
            information in the ATOM records of the PDB entry and 
            the sequence database entry given on DBREF. Please note 
            that these records were designed to identify differences
            and not errors.  

            COLUMNS        DATA TYPE       FIELD          DEFINITION                          
            ----------------------------------------------------------------------------------
             1 -  6        Record name     "SEQADV"                                           
             8 - 11        IDcode          idCode         ID code of this entry.              
            13 - 15        Residue name    res_name        Name of the PDB residue in conflict.
            17             Character       chainID        PDB chain identifier.               
            19 - 22        Integer         seqNum         PDB sequence number.                
            23             AChar           iCode          PDB insertion code.                 
            25 - 28        LString         database       Sequence database name.             
            30 - 38        LString         dbIdCode       Sequence database accession     
                                                          number.                             
            40 - 42        Residue name    dbRes          Sequence database residue name.     
            44 - 48        Integer         dbSeq          Sequence database sequence number.  
            50 - 70        LString         conflict       Conflict comment. 
            
            List of comments:
                    Cloning artifact, Expression tag, Conflict, Engineered Variant,
                    Insertion, Deletion, Microheterogeneity, Chromophore,
                    See remark 999

        """
        kw = self.type_to_remark['seqAdv']
        new = {}
        for l in data:
            mod_desc = l[49:].strip()
            if "deletion" in mod_desc.lower():
                if ignore_deletion:
                    continue
            pdb_id = l[7:11].strip()
            res_name = l[12:15].strip()
            ch_id = l[16] #.strip()
            try:
                res_num = int(l[18:22].strip())
            except:
                res_num = None
            iCode = l[22].strip()
            std_res_name = l[39:42].strip()
            #if not ch_id in new:
            #    new[ch_id] = []
            chain = new.setdefault(ch_id, {})
            if not res_num in chain:
                chain[res_num] = {'res_name': res_name,
                'res_num' : res_num, 'iCode' : iCode,
                'std_res_name' : std_res_name, 'mod_desc' : mod_desc}
            else:
                if not type(chain[res_num]['std_res_name']) == list:
                    chain[res_num]['std_res_name'] = [chain[res_num]['std_res_name']] 
                chain[res_num]['std_res_name'].append(std_res_name) 
        return new

    def _parse_symmetry(self, data):
        # key = 'REMARK 290'
        # terminated by 'REMARK 290 REMARK: NULL'
        return data

    def _parse_biounit(self, data):
        """ Test with Stout's protease structure!"""
        return data
        data = new

        """
        REMARK 350 BIOMOLECULE: ?
        REMARK 350 APPLY THE FOLLOWING TO CHAINS: ?, ?...
        REMARK 350   BIOMT1   N  N.NNNNNN  N.NNNNNN  N.NNNNNN        N.NNNNN
        REMARK 350   BIOMT2   N  N.NNNNNN  N.NNNNNN  N.NNNNNN        N.NNNNN
        REMARK 350   BIOMT3   N  N.NNNNNN  N.NNNNNN  N.NNNNNN        N.NNNNN
        """

    def _parse_het_list(self, data):
        """ 
        [ SOURCE:http://www.wwpdb.org/documentation/format23/sect4.html ]

        COLUMNS     DATA TYPE     FIELD         DEFINITION
        ------------------------------------------------------
         1 -  6     Record name   "HET      "
         8 - 10     LString(3)    hetID         Het identifier, right-justified.
        13          Character     ChainID       Chain identifier.
        14 - 17     Integer       seqNum        Sequence number.
        18          AChar         iCode         Insertion code.
        21 - 25     Integer       numHetAtoms   Number of HETATM records for the
                                                group present in the entry.
        31 - 70     String        text          Text describing Het group.
        """
        new = {}
        kw = self.type_to_remark['hetList']
        for l in data:
            het_id = l[7:10].strip()
            chain_id = l[12].strip()
            seq_num = l[13:17].strip()
            iCode = l[17]
            num_het_atoms = l[20:25].strip()
            text = l[30:70].strip()
            if not chain_id in list(new.keys()):
                new[chain_id] = []
            new[chain_id].append( (het_id,seq_num) )
        return new


    def init_vars(self):
        """ """
        self._kwslicing =[0,6]
        self.kw = { 'ATOM  ':None, 
                    'HETATM':None,
                    'MODEL ':None,
                    'ENDMDL':None,
                    'TER   ':None,
                   #'MASTER ':None, # XXX Neither is supported! 
                   #'END   ':None,  # XXX
                    'CONECT':None, 
                  }

        self.coord_kw = ('ATOM  ', 'HETATM')

        self.model_set = {}

        # list of specific altmodes per residue
        self.alt_res_mode = {}

        # list of altresidues found per model
        self.alt_residues = {}

        # default altLocation mode
        self.alt_mode = 'A'

        self.conect = []
        self.header = []
        # TODO implement the CONECT to find if multi-chains are linked
        # TODO also look at LINK kw that provides *EXPLICIT* links between chains (see 1o7d)
        """ 
        [ source: http://www.wwpdb.org/documentation/format33/sect10.html#CONECT ]
        COLUMNS       DATA  TYPE      FIELD        DEFINITION
        -------------------------------------------------------------------------
         1 -  6        Record name    "CONECT"
         7 - 11        Integer        serial       Atom  serial number
        12 - 16        Integer        serial       Serial number of bonded atom  _
        17 - 21        Integer        serial       Serial number of bonded atom   |
        22 - 26        Integer        serial       Serial number of bonded atom   +-- optional
        27 - 31        Integer        serial       Serial number of bonded atom  _| 
        """

    def parse_lines(self):
        """ does the dirty work parsing header and models """
        inside = False
        curr_model = []
        i,j = self._kwslicing
        # model state of the molecule
        self.current_model = 0
        model_id = 0
        for idx, raw in enumerate(self.text):
            l = raw[i:j]
            if not l in self.kw:
                self.header.append(idx)
                continue
            if l == 'ENDMDL':
                inside = False
                curr_model.append(idx)
                self.model_set[model_id] = curr_model
                model_id += 1
                curr_model = []
                #print "FLUSHING", curr_model, model_id
            elif l == 'MODEL ':
                inside = True
                curr_model = [idx]
            elif l == 'CONECT':
                pass
                #bonds = self._parseConect(raw)
                #self.conect.append(bonds)
            else:
                curr_model.append(idx)
        if len(curr_model):
            self.model_set[model_id] = curr_model


    def set_state(self, model=None, alt_mode=None, alt_res_mode={}):
        """ define the state of the multistructure"""
        if model:
            self.current_model = model
        self.set_alt_mode(alt_mode)
        self.set_res_alt_mode(alt_res_mode)
        #self.generate_obmol()


    def set_alt_mode(self, alt_mode='a'):
        """ define the default alternate conformation mode for
            residues that do not have a specific setting in self.alt_res_mode
        """
        if not alt_mode:
            return
        self.alt_mode = alt_mode.upper()

    def set_res_alt_mode(self, alt_res_mode={}):
        """ set the residue-specific alternate residue mode"""
        if alt_res_mode == {}:
            self.alt_residue_mode = {}
            return
        for k,v in alt_res_mode:
            self.alt_res_mode[k] = v

    def generate_obmol(self, raw):
        """ create the OBMolecule from the current model"""
        self.mol = ob.OBMol()
        conv = ob.OBConversion()
        conv.SetInFormat('pdb')
        conv.ReadString(self.mol, raw)
        return self.mol

    def get_header(self):
        """ return PDB header information """
        return [self.text[i] for i in self.header]


    def parse_atoms(self):
        """ process all atoms found in each model"""
        #self.graph = []
        # XXX
        kw = self.coord_kw # ('ATOM   ', 'HETATM ')
        i,j = self._kwslicing
        for mIdx, model in list(self.model_set.items()):
            self.current_model = mIdx
            self.alt_residues[mIdx] = {}
            structure = {}
            # scan all atom lines in the model
            for line_idx in model:
                raw = self.text[line_idx]
                if raw[i:j] in kw:
                    self.add_atom_to_struct(structure, line_idx)
            # cleanup the model
            clean_structure = self.clean_structure(structure)
            # store the model
            self.model_set[mIdx] = clean_structure
        # reset model to the begin
        self.current_model = 0
        
        # DEBUG
        if 0:
            for ch, re in list(self.model_set[mIdx].items()):
                print("CH[%s]" % ch)
                for nu, inf in list(re.items()):
                    for k,v in list(inf.items()):
                        print("\t", k, v)
                


    @property
    def model_count(self):
        """ return the number of models found in the structure"""
        return len(self.model_set)

    @property
    def get_alt_residues(self):
        # TODO it was get_alt_residues
        """ return the number of models found in the structure"""
        return sorted(self.alt_residues[self.current_model].keys())

    @ property
    def chains(self):
        return sorted(self.model_set[self.current_model].keys())

    def add_atom_to_struct(self, structure, line_idx):
        """ take care of creating/adding atoms, book-keeping..."""
        chain_info, res_info, at_info = self.get_atom_info(line_idx)
        res_name = res_info['name']
        res_num = res_info['num']
        chain = chain_info['name']
        #resKey = "%s%s" % (res_name, res_num)
        # it is possible that an alternate location defines a different residue!
        # therefore the key is going to be the sequence number and not the NameNumber
        resKey = "%d" % (res_num)
        residues = structure.setdefault(chain, {})
        atoms = residues.setdefault(resKey, [] )
        atoms.append(at_info)


    def clean_structure(self, structure):
        """process each residue for alt states"""
        for chain, residue in list(structure.items()):
            for res, atoms  in list(residue.items()):
                new_res_atoms = self.compact_atoms(atoms)
                structure[chain][res] = new_res_atoms
        return structure


    def compact_atoms(self, atoms):
        """ generate the structure to create requested alt conformations
        """
        alt_residues = self.alt_residues[self.current_model]
        common = []
        alt_atoms = {}
        for a in atoms:
            if a['is_alt']:
                alt_label = a['alt']
                # register residue in list of alt residues 
                alt_list = alt_residues.setdefault(a['res_num'], [])
                if not alt_label in alt_list:
                    alt_list.append(alt_label)
                # register atoms to the specific alt group
                alt_group = alt_atoms.setdefault(alt_label, [])
                alt_group.append(a)
            else:
                common.append(a)
        return {'common': common, 'alternate': alt_atoms}


    def get_structure(self, chains=[]):
        """ generate OBMol molecule structure (optionally, 
            containing only selected chains
        """
        raw = "".join( self.get_struct_state(chains=chains))
        #open('DEBUG_XX.pdb','w').write(raw) 
        return self.generate_obmol(raw)
        
        

    def get_struct_state(self, chains=[]): #, alt_mode='A', model='all', alt_res_mode=[]):
        """ create a representation o the structure using 
        
            the current model (self.model), alternate locations
            per residue (self.alt_res_mode) and default (self.alt_mode)
        """
        out = []
        curr_model = self.model_set[self.current_model]
        for model_chain, residues in list(curr_model.items()):
            if chains and not model_chain in chains:
                continue
            for res_num in sorted(residues.keys()):
                res_atoms = residues[res_num]
                res_state = self.get_res_state(res_num, res_atoms)
                out += res_state
        return [self.text[i] for i in out]

    def sorted_residues(self, residues):
        """ return a tuple with """
        out = []
        keys = list(residues.keys())
        name_num = [ (x[0:3], int(x[3:])) for x in keys ]
        name_num.sort(key=itemgetter(1))
        for name, num in name_num:
            k = "%s%d" % (name, num)
            out.append( (k, residues[k]) )
        return out
        
    def get_res_state(self, res_num, res_atoms):
        """ extract residue atoms in the specified alt_modes (default or specific)"""
        # common atoms with no alternate conformations
        out = [ x['line_idx'] for x in res_atoms['common'] ]

        alt_mode = self.alt_res_mode.setdefault(res_num, self.alt_mode)
        if alt_mode in res_atoms['alternate']:
            for a in res_atoms['alternate'][alt_mode]:
                    out.append(a['line_idx'])
        return out


    def get_atom_info(self, line_idx):
        """ extract information from the ATOM/HETATM line"""
        # atom information
        s = self.text[line_idx]
        atName = s[12:16].strip()
        atNum = int(s[6:11])
        alt = s[16].strip()
        is_alt = bool(len(alt))
        #occupancy = float(s[54:60].strip())
        #temp = float(s[60:66].strip())
        #element = s[76:78]
        #coord = map(float, ( s[30:38], s[38:46], s[46:54]) )

        # remove the alt location label
        self.text[line_idx] = s[:16] + " "+ s[17:]

        # chain
        chain = s[21]
        #segment = s[72:76]
        chain_info = { 'name': chain}
        # residue
        res_num = int(s[22:26].strip())
        res_name = s[17:20].strip()
        res_id = "%s:%s%d" % (chain, res_name, res_num)
        res_info = {'name': res_name, 'num': res_num, 'res_id':res_id}
        # atom
        at_info = { 'is_alt' : is_alt, 'atom_name': atName, 
            'atom_num': atNum, 
            'line_idx' : line_idx, 'res_id': res_id, 'alt':alt, 'res_num':res_num,
            #'segment':segment,
            #'element':element,
            #'occupancy': occupancy, 'temp':temp, 'coord':coord,
            }
        return chain_info, res_info, at_info

                
        AaaToA = { 'gly' : 'G', 'ala': 'A', 'val': 'V', 'leu': 'L', 'ile':'I', 'met':'M',
            'phe': 'F', 'trp' : 'W', 'pro': 'P', 'ser':'S', 'thr':'T', 'cys':'C', 'tyr':'Y',
            'asn':'N', 'gln': 'Q', 'asp':'D', 'glu':'E', 'lys':'K', 'arg':'R', 'his':'H'}



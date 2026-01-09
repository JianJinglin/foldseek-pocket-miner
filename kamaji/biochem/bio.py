########################################################################
#
# Date: 2016 Authors: Stefano Forli
#
#    forli@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Stefano Forli and TSRI 2016
#
#########################################################################


# strictly PDB?

aa_3l_to_1l = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


aa_1l_to_3l = {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY',
                        'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M': 'MET',
                        'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER',
                        'R': 'ARG', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR'}

aa_atom_count = { 'ALA': 6,  # resName : heavy_atoms_count
      'ARG': 12, 'ASN':9, 'ASP':9, 'ASX':9,
      'CYS':7, 'GLU':10, 'GLN':10, 'GLX':10, 'GLY':5,
      'HIS':11, 'ILE':9, 'LEU':9, 'LYS':10, 'MET':9,
      'PHE':12, 'PRO':8, 'SER':7, 'THR':8, 'TRP':15,
      'TYR':13, 'VAL':8 }

aa_bond_count = {'HIS': 10, 'ILE': 7, 'LEU': 7, 'LYS': 8,
        'MET': 7, 'PHE': 11, 'PRO': 7,
        'SER': 5, 'THR': 6, 'TRP': 15,
        'ALA': 4, 'TYR': 12, 'VAL': 6,
        'ARG': 10, 'ASN': 7, 'ASP': 7,
        'CYS': 5, 'GLN': 8, 'GLU': 8,
        'GLY': 3}

aa_data_bb = {
    # Adaptation from original from DayLight:
    #     http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
    # extra O/N at the end have been removed, because troublesome.
    #'pro' : '[$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]',
    #'generic' : '[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])',
    #'gly' : '[$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N])]',
    #
    # Also, Lys and Arg patterns have been modified to have an extra generic N, otherwise
    # OpenBabel is not able to match them.
    # modified from origianl Daylight to include
    # 'generic' : '[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[$[OX2H,OX1-,N]]',
    'generic' : '[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])',
    'generic_patternIndices' : { 'N':0, 'CA':1, 'C':2, 'O':3 },
    # do not match Pro and Gly
    'pro' : '[$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])',
    # 'gly' : '[$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1]))]',
    #  fixed from daylight (the pattern matches only a single atom)
    'gly' : '[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]'
        }

backbone = { 'labels' : [ ['N','CA'], ['C', 'O'] ]}

proline = [ 'N', 'CA', 'CB', 'CG', 'CD', 'C', 'O' ]

aa_data_ss = {
    'ala': { 'pattern' : '[CH3X4]',
            'labels' : ['CB'], },

    'arg': { 'pattern' : '[CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0,N])[NH2X3]',
            # Hits acid and conjugate base. MODIFIED from original implementation
            # adding an extra generic 'N' type for the head, otherwise OB will
            # miss it
             'labels' : [ 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2' ],
             },

    'asn': { 'pattern' : '[CH2X4][CX3](=[OX1])[NX3H2]', #
            # Also hits Gln side chain when used alone
             'labels' : [ 'CB', 'CG', 'OD1','ND2' ],
            },

    'gln': { 'pattern' : '[CH2X4][CH2X4][CX3](=[OX1])[NX3H2]', #
            # Also hits Gln side chain when used alone
             'labels' : [ 'CB', 'CG', 'CD', 'OE1','NE2' ],
            },

    'asp': { 'pattern' : '[CH2X4][CX3](=[OX1])[OH0-,OH]',
            # Aspartate (or Aspartic acid) side chain. Hits acid and conjugate base.
            # Also hits Glu side chain when used alone.
             'labels' : [ 'CB', 'CG', 'OD1','OD2' ],
             },

    'cys' : { 'pattern': '[CH2X4][SX2H,SX1H0-]', # Cysteine side chain. Hits acid and conjugate base
             'labels' : [ 'CB', 'SG' ],
             },

    'glu' : { 'pattern' : '[CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]', # Hits acid and conjugate base
              'labels' : [ 'CB', 'CG', 'CD', 'OE1', 'OE2'] ,
              },

    'his' : { 'pattern' : ('[CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:'
                           '[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1'),
            #Hits acid & conjugate base for either Nitrogen. Note that the Ns can be either
            # ([(Cationic 3-connected with one H) or (Neutral 2-connected without any Hs)]
            # where there is a second-neighbor who is [3-connected with one H]) or (3-connected with one H).
            'labels': ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
            },

    'ile' : { 'pattern' : '[CHX4]([CH3X4])[CH2X4][CH3X4]',
            'labels' : ['CB', 'CG1', 'CG2', 'CD1'] ,
            },


    'leu' : { 'pattern' : '[CH2X4][CHX4]([CH3X4])[CH3X4]',
             'labels' : ['CB', 'CG', 'CD1', 'CD2' ]
             },

    'lys' : { 'pattern' : '[CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0,N]',
            # MODIFIED from original implementation
            # adding an extra generic 'N' type for the head, otherwise OB will
            # miss it
             'labels' : [ 'CB', 'CG', 'CD', 'CE', 'NZ' ],
             },

    'met' : { 'pattern': '[CH2X4][CH2X4][SX2][CH3X4]',
             'labels' : [ 'CB', 'CG', 'SD', 'CE' ],
             },


    'phe' :  { 'pattern' : '[CH2X4][cX3](1[cX3H][cX3H][cX3H][cX3H][cX3H]1)',
              'labels' : [ 'CB', 'CG', 'CD1', 'CD2', 'CZ', 'CE1', 'CE2' ] ,
              },

    'ser' :  {'pattern': '[CH2X4][OX2H]',
                'labels' : [ 'CB', 'OG' ],
                },

    'thr' : {'pattern' : '[CHX4]([CH3X4])[OX2H]',
             'labels' : [ 'CB', 'CG2', 'OG1' ],
             },

    'trp' : { 'pattern' : '[CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12',
             'labels' : [ 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2' ],
             },

    'tyr' : {'pattern': '[CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1', # Acid and conjugate base
            'labels' : [ 'CB', 'CG', 'CD1', 'CD2', 'CZ', 'OH', 'CE1', 'CE2'],
            },

    'val' : { 'pattern' : '[CHX4]([CH3X4])[CH3X4]',
            'labels' : [ 'CB', 'CG1', 'CG2'] ,
            },
        }

"""
mask format: gets applied to the smart pattern:
    0 : ignore
   -1 : delete
    n : change to element n
LIMITATIONS:
  - working with mol copies, it doesn't handle nicely atoms
    deleted at the interface with other residues
smarts have been created from here
https://cactus.nci.nih.gov/chemical/structure/[SMILES HERE]/file?%20format=smarts
"""

# TODO possibly obsolete
resModifications = {
    'MSE': { 'from': { 'smi':'C[Se]CC[C@H](N)C(O)=O',
                       'smarts': '[#6]-[Se]-[#6]-[#6]-[#6@H](-[#7])-[#6](~[#8])', #=[#8]',
                       'mask': (0, 16, 0, 0, 0, 0, 0, 0),
                       'name':'selenomethionine'},
             'to' :  {'smi': 'CSCC[C@H](N)C(O)=O',
                       'smarts' : '[#6]-[S]-[#6]-[#6]-[#6@H](-[#7])-[#6](-[#8])=[#8]',
                       'name' : 'MET'},
           },

    'SEP': { 'from': {'smi':'N[C@@H](COP(O)(O)=O)C(O)=O',
                    'name': 'phosphoserine'},
             'to' :  {'smi': 'N[C@@H](CO)C(O)=O',
                      'name' : 'SER'},
                      },

    'CAS': { 'from': {'smi':'C[As](C)SC[C@H](N)C(O)=O',
                      'smarts': '[#6]-[As](-[#6])-[#16]-[#6]-[#6@H](-[#7])-[#6](~[#8])', #=[#8]',
                       'mask': (-1, -1, -1, 0, 0, 0, 0, 0,0),
                      'name': 'S-(dimethylarsenic)cysteine'},
             'to' :  {'smi': 'N[C@@H](CS)C(O)=O',
                      'smarts': '[#16]-[#6]-[#6@H](-[#7])-[#6](-[#8])=[#8]',
                      'name' :'CYS'}, },

    'TPO' : {'from': {'smi' :'C[C@@H](OP(O)(O)=O)[C@H](N)C(O)=O',
                      'name' : 'phosphothreonine'},
             'to': {'smi': 'C[C@@H](O)[C@H](N)C(O)=O',
                      'name' : 'THR'},
            },
    # TODO ADD post-translational mods (ABA, methylation, kcx, all l-peptide linking, etc.)
}


dnabasesAtomCount = { 'DA' : 22,  # resName : heavy_atoms_count
                'DC' : 20, 'DG': 23, 'DT': 21 }

dnabasesBondCount = { 'DA': 23, 'DC': 20, 'DG': 24, 'DT': 21 }

rnabasesAtomCount = {'A':23, # resName : heavy_atoms_count
        'C':21, 'G':24, 'U':21}

rnabasesBondCount = { 'A': 24, 'C': 21, 'G': 25, 'U': 21}

waterNames = [ 'HOH', 'WAT', 'DOD' ]

cofactors = {
            'HEM': { 'name': 'Protoporphyrin IX containing Fe (Heme)',
                    'smi':  ('CC1=C(CCC(O)=O)C2=Cc3c(CCC(O)=O)c(C)c4C='
                               'C5C(C)=C(C=C)C6=[N]5[Fe]5([N]2=C1C=c1c'
                               '(C=C)c(C)c(=C6)n51)n34'),
                    'setup'  : ['metal'] }, # defines the kind of setup required for docking?

###             'HEC': { 'name': 'Heme C',
###                 'smi':  ('C\C=C1/C(=C2C=C3N4C(=Cc5n6c(C=C7N8C(=C(C)\C7=C/C)'
###                         'C=C1N2[Fe@@]468)c(C)c5CCC(O)=O)C(=C3C)CCC(O)=O)C'),
###                     'setup'  : ['metal'] },
###

            'ADP' : {'name' : 'Adenosine di-phosphate (ADP)',
                     'smi':  ('Nc1ncnc2n(cnc12)[C@@H]1O[C@H]'
                              '(CO[P@@](O)(=O)OP(O)(O)=O)[C@@H]'
                              '(O)[C@H]1O'),
                    'setup' : []},

            'ATP' : {'name' : 'Adenosine tri-phosphate (ATP)',
                     'smi': ('Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO[P@](O)'
                             '(=O)O[P@@](O)(=O)OP(O)(O)=O)[C@@H](O)'
                             '[C@H]1O'),
                     'setup' : []},

            'NAD' : {'name': 'Nicotinamide adenine dinucleotide (NADP)',
                     'smi':  ('NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H]'
                              '(CO[P@]([O-])(=O)O[P@](O)(=O)OC'
                              '[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2'
                              'cnc3c(N)ncnc23)[C@@H](O)[C@H]1O'),
                    'setup' : []},

            'NDP' : {'name': 'Dihydro-nicotinamide-adenine-dinucleotide (NADP-ph)',
                     'smi':  ('NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](CO[P@@]'
                              '(O)(=O)O[P@](O)(=O)OC[C@H]2O[C@H]([C@H]'
                              '(OP(O)(O)=O)[C@@H]2O)n2cnc3c(N)ncnc23)'
                              '[C@@H](O)[C@H]1O'),
                    'setup' : []},

            'NAI' : {'name': 'Nicotininamide adenin dinucleotide (NADPH+)',
                    'smi' :  ('NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H]'
                              '(CO[P@@](O)(=O)O[P@](O)(=O)OC'
                              '[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2'
                              'cnc3c(N)ncnc23)[C@@H](O)[C@H]1O'),
                    'setup' : []},

            'GDP' : { 'name' : "Guanosine-5'-diphosphate",
                      'smi' : ( 'Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H]'
                            '(CO[P@@](O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O'),
                    'setup' : []},

            'GTP' : { 'name' : "Guanosine-5'-triphosphate",
                     'smi' : ('Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O'
                              '[C@H](CO[P@@](O)(=O)O[P@@](O)(=O)'
                              'OP(O)(O)=O)[C@@H](O)[C@H]1O'),
                    'setup' : []},

            'PLP' : { 'name' : "Pyridoxal-5'-phosphate (vitamin B6 phosphate)",
                      'smi'   : 'Cc1ncc(COP(O)(O)=O)c(C=O)c1O',
                      'setup' : []},

            'FAD' : { 'name': "Flavin-adenine dinucleotide",
                      'smi':  ('Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C[C@H]'
                               '(O)[C@H](O)[C@H](O)CO[P@](O)(=O)O[P@@]'
                               '(O)(=O)OC[C@H]3O[C@H]([C@H](O)[C@@H]3O)'
                               'n3cnc4c(N)ncnc34)c2cc1C'),
                    'setup' : []},

            'FMN' : { 'name': "Flavin-adenine mononucleotide",
                      'smi':  ('Cc1cc2nc3c(nc(=O)[nH]c3=O)n'
                                '(C[C@H](O)[C@H](O)[C@H](O)'
                                'COP(O)(O)=O)c2cc1C'),
                    'setup' : []},

            'SF4' : { 'name': 'Iron/sulfur cluster',
                      'smi' : ('[S@@]12[Fe]3[S@]4[Fe]5[S@@'
                                    ']3[Fe]1[S@]5[Fe]24'),
                      'setup': ['metal'] },

            'CR2' : { 'name': 'GPF chromophore',
                      'smi' : ('NCC1=N\C(=C/c2ccc(O)cc2)C(=O)N1CC(O)=O'),
                      'setup': [''] },

                            }

additives = {  # crystallography additives/buffers
        'SO4' : { 'name' : 'sulfate' , 'smi' : '[O-]S([O-])(=O)=O', 'role':'buffer'},
        'GOL' : { 'name' : 'glycerol' , 'smi': 'C(C(CO)O)O', 'role':'additive'},
        'EDO' : { 'name' : 'ethylen glycol' , 'smi': 'OCCO', 'role':'additive'},
            #'NAG' : { 'name'
        'SDS' : { 'name' : 'dodecyl sulfate' , 'smi': 'CCCCCCCCCCCCOS(O)(=O)=O'} ,
        'TRS' : { 'name' : 'tris-buffer' , 'smi': '[NH3+]C(CO)(CO)CO', 'role':'buffer'},
        'LDA' : { 'name' : 'Lauryl dimethylamine-N-oxide' , 'smi': 'CCCCCCCCCCCC[N+](C)(C)[O-]'},
        'DMS' : { 'name' : 'DMSO' , 'smi': 'CS(C)=O', 'role':'additive'},
            #'DM2' : { 'name' : 'DMSO (low quality)' , 'smi': 'CS(C)O'},
        'PO4' : { 'name' : 'phosphate', 'smi': '[O-]P([O-])([O-])=O','role':'buffer'},
        'ACT' : { 'name' : 'acetate' , 'smi': 'CC([O-])=O','role':'buffer'},
        'EPE' : { 'name' : 'HEPES' , 'smi': 'OCCN1CCN(CCS(O)(=O)=O)CC1', 'role':'buffer'},
        'NH4' : { 'name' : 'ammonium' , 'smi': '[NH4+]','role':'buffer'},
        'PEG' : { 'name' : 'PEG (poly[Ethylen-glycole])' , 'smi': 'OCCOCCO', 'role':'additive'},
        'EOH' : { 'name' : 'ethanol' , 'smi': 'CCO', 'role':'additive'},
        'BME' : { 'name' : 'Beta-mercapto-ethanol' , 'smi': 'OCCS', 'role':'reducing agent'},
        'FMT' : { 'name' : 'formic acid' , 'smi': 'OC=O', 'role':'buffer'},
        'PYR' : { 'name' : 'pyruvic acid' , 'smi': 'CC(=O)C(O)=O', 'role':'buffer'},
        # too many of these, a SMARTS should be better...
        # { 'name' : 'tartaric acid', 'smi': 'O[C@H]([C@@H](O)C(O)=O)C(O)=O'},
        'MES' : { 'name' : '2-morpholin-4-ium-4-ylethanesulfonate' , 'smi': '[O-]S(=O)(=O)CC[NH+]1CCOCC1' , 'role':'buffer'},
        'BCT' : { 'name' : 'bicarbonate' , 'smi': 'OC([O-])=O', 'role':'buffer'},
        'NHE' : { 'name' : '2-[N-cyclohexylamino]ethane sulfonic acid', 'smi': 'O[S](=O)(=O)CCNC1CCCCC1', 'role':'buffer' },
        #'PEG' : { 'name' : 'PEG (poly[Ethylen-glycole])' , 'smarts': '[O][C][C][O]', 'role':'additive'},
            }


residuenormalization =  { 'MSE': { 'from': { 'smi':'C[Se]CC[C@H](N)C(O)=O',
                                              'name':'selenomethionine'},
                                     'to' :  {'smi': 'CSCC[C@H](N)C(O)=O',
                                              'name' : 'MET'}, },

                            'SEP': { 'from': {'smi':'N[C@@H](COP(O)(O)=O)C(O)=O',
                                            'name': 'phosphoserine'},
                                     'to' :  {'smi': 'N[C@@H](CO)C(O)=O',
                                              'name' : 'SER'}, },

                            'CAS': { 'from': {'smi':'C[As](C)SC[C@H](N)C(O)=O',
                                              'name': 'S-(dimethylarsenic)cysteine'},
                                     'to' :  {'smi': 'N[C@@H](CS)C(O)=O',
                                              'name' :'CYS'}, },

                            'TPO' : {'from': {'smi' :'C[C@@H](OP(O)(O)=O)[C@H](N)C(O)=O',
                                              'name' : 'phosphothreonine'},
                                     'to': {'smi': 'C[C@@H](O)[C@H](N)C(O)=O',
                                              'name' : 'THR'}, },
                          }





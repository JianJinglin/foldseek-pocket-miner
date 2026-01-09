
# Description of properties of potentially coordinating metals
#
#
#
#                'Fe' : { 
#                            0 : { 'coordNum'   :   2,     # how many atoms could coordinate max (including res,lig,wat...)
#                               'dist'       :   2.5,    # coordinating distance (add +0.3 for S)
#                               'text'       :   ('Fe2+ as found in heme'),     # text description
#                               'minRes'     :   4,  # min.count of residues (i.e., non-ligand, non-water) 
#                                                    # to be considered coordinating and not structural
#                                #'charge' :  0.75 * 2,
#                                },
#                            1: { ... }
#
#                        },
#
                            



# IRON AND MIDDLE-WATER
# check this paper: http://ac.els-cdn.com/S0022283609002320/1-s2.0-S0022283609002320-main.pdf?_tid=e12deee8-1ac7-11e5-bf8d-00000aab0f6b&acdnat=1435188328_3137ab1e37691a8505595a1fd6694ab3
metals = [ 'Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 
                      'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
                      'Ni', 'Cu', 'Zn', 'Ga', 'Ge',  # 'As', 'Se', 
                      'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 
                      'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
                      'Sb', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 
                      'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 
                      'Pb', 'Bi', 'Po', 'Fr', 'Ra', 'Rf', 
                      'Db', 'Sg', 'Bh', 'Hs', 'La', 'Ce', 'Pr', 
                      'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
                      'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th',
                      'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 
                      'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr' ]
 

nonmetals = metalloids = {'Si':14, 'P':15, 'S':16, 'Se':34, 'As':33}

coordinatingMetals = {
                'Zn' : { 'dist': 2.5,  'coordinatingElements': ['O', 'S', 'N' ],
                        'atomicNum' : 30,
                         'contacts': { 3 : { 'text' : ('Zn2+ tetrahedral coordination'),
                                              'catalytic': True , },
                                       4 : { 'text'      :   ('Zn2+ with structural role'),
                                               'catalytic': False },
                                       2 : { 'text' : ('Zn2+ possible artifact (i.e., crystallization additive)'),
                                              'catalytic': False }, },
                                                        },
                'Fe' : { 'dist' : 2.5,  'coordinatingElements': ['O', 'S', 'N' ],
                        'atomicNum' : 26,
                        'contacts': {  5 : { 'text' : ('Fe2+ coordinating as found in heme'),
                                           'catalytic': True },
                                       6 : { 'text'       :   ('Fe2+ non-coordinating (heme)'),
                                          'catalytic': False },
                                       4 : { 'text'       :   ('Fe2+ non-coordinating inorganic (iron cluster, detached heme)'),
                                           'catalytic': False },
                                       3 : { 'text'       :   ('Fe2+ coordinating inorganic (iron cluster)'),
                                            'catalytic': True }, },
                                                        },
                'Mg' : { 'dist': 2.4,  'coordinatingElements': ['O', 'S', 'N' ], 
                        'atomicNum' : 12,
                         'contacts' : { 2 : { 'text' : ('Mg+2 mono/bidentate coordinated'),
                                          'catalytic': True },
                                       3 : { 'text'       :   ('Mg+2 mono/bidentate coordinated'),
                                          'catalytic': True }, },
                                          },
                'Mn' : { 'dist': 2.4, 'coordinatingElements': ['O', 'S', 'N' ], 
                        'atomicNum' : 25,
                         'contacts' : { 4 : { 'text' : ('Mn+2 coordinated'),
                                        'catalytic': True },
                                        2 : { 'text'       :   ('Mn+2 coordinated'),
                                            'catalytic': True }, }
                                        },
                        }


catalyticMetalsOld = {
                'Zn' : { 0: {'coordNum'  :   3,  # receptor coordination 'points'
                            'dist'      :   2.5,
                            'text'      :   ('Zn2+ tetrahedral coordination'),
                            'coordinating': True },
                        1: {'coordNum'   :   4,       # receptor coordination points
                            'dist'      :   2.5,
                            'text'      :   ('Zn2+ with structural role'),
                            'coordinating': False },
                        2: {'coordNum'     :   2,
                            'dist'      :   2.5,
                            'text'      :   ('Zn2+ possible artifact (i.e., crystallization role)'),
                            'coordinating': False }, },
                

                'Fe' : { 0 : { 'coordNum'   :   5,     # how many atoms could coordinate max
                               'dist'       :   2.5,    # coordinating distance
                               'text'       :   ('Fe2+ coordinating as found in heme'),
                               'coordinating': True },
                         1 : {  'coordNum'   :   6,     # how many atoms could coordinate max
                               'dist'       :   2.5,    # coordinating distance +0.3 for S
                               'text'       :   ('Fe2+ non-coordinating as found in heme'),
                               'coordinating': False },
                         2 : {  'coordNum'   :   4,     # how many atoms could coordinate max
                               'dist'       :   2.5,    # coordinating distance +0.3 for S
                               'text'       :   ('Fe2+ non-coordinating inorganic (iron cluster)'),
                               'coordinating': False },
                         3 : {  'coordNum'   :   3,     # how many atoms could coordinate max
                               'dist'       :   2.5,    # coordinating distance +0.3 for S
                               'text'       :   ('Fe2+ coordinating inorganic (iron cluster)'),
                               'coordinating': True }, },

                'Mg' : { 0 : { 'coordNum'   :   2,
                               'dist'       :   2.3,
                               'text'       :   ('Mg+2 mono/bidentate coordinated'),
                                'coordinating': True },
                         1 : { 'coordNum'   :   3,
                               'dist'       :   2.3,
                               'text'       :   ('Mg+2 mono/bidentate coordinated'),
                                'coordinating': True }, },

                'Mn' : { 0 :{ 'coordNum'   :   4,
                               'dist'       :   2.4,
                               'text'       :   ('Mn+2 coordinated'),
                                'coordinating': True },

                         1 :{ 'coordNum'   :   2,
                               'dist'       :   2.4,
                               'text'       :   ('Mn+2 coordinated'),
                                'coordinating': True },
 
                            }
                }
                    
                    
                    
                # Ca: 4:1b8l 6: catalytic?
                    
                #    , 'Co', 'Ni', 'Cu'}

salts = [ 'Na', 'Cl', 'Ca', 'K', 'Li', 'Cl', 'N' ]

simpleorganics = {
        # glycane pattern ( O_inring *RINGBOND* C_inring *RINGBOND*
        #                   nonRing_C_ *NONRINGBOND*-NonRingOxygen )
        #'glycosilation':'[#8r5,#8r6]@[#6r5,#6r6]~[C!r]~[#8!r]',
        # TODO fix to capture also HSY
        # FIXME this does not capture also FUC!!!
        # TODO add two patterns for rings?
        'glycosilation':'[#8r5,#8r6]@[#6r5,#6r6]~[C!r]~[#8!r]',
        'lipids' : [ # <- list because order is important
              # phospholipid (phosphatidyl-choline group)
              #{'phosphatidic acid': 'C(OC=O)(COC=O)COP(~O)(~O)~O'},
              # carboxylic-lipid pattern with at least 4 carbons attached to at least one hydrogen
              {'phospholipid': 'C(OC~O)(COC~O)COP(~O)(~O)~O'}, # NO BOND ORDER
              #{'lipid acid': '[C;!R;h1]~[C;!R;h1]~[C;!R;h1]~[C;!R;h1]~[CX3](=[OX1])O'},
              {'fatty acid': '[C;!R;H2]~[C;!R;H2]~[C;!R;H2]~[C;!R;H2]~[CX3](~[OX1])O'},
              # {'fatty acid SMI': 'CCCCC(=O)O'},
              # aliphatic chain with 6-any carbons attached to an oxygen
              #{'lipid (generic)': '[Ch1;!R]~[Ch1;!R]~[Ch1;!R]~[Ch1;!R]~[CX3][#8].[!r]'},
              #{'lipid (generic)': '([Ch1;!R]~[Ch1;!R]~[Ch1;!R]~[Ch1;!R]~[CX3][#8]),([*!r])'},
              #{ 'lipid (generic)' :  '[C]~[C]~[C]~[C]~[C]~[#8]'},
              ]
        }

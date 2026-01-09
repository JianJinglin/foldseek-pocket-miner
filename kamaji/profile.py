
# from Raccoon.utils.debugtools import DebugObj


# properties
# .nucleic
# .protein
# .dna
# .rna
# .modified
# .incomplete
# .chains
# .water
# .ions
# .cofactor
# .additive
# .modifier
# .ligand
# .salt
# .metal
# .catalytic
#

# for c in container.chains:
#   str(c) == True?

_required = [ 'chain', 'name', 'num', 'id' ]
#_required = [ 'chain' ]

from pprint import pprint as pp

class StructureProfile(object):
    """ This class is used to manage the information collected by Kamaji
        about a target structure.
        >>> k = StructureProfile()
        >>> k.add_residue( OBResidue.GetIdx(), 'protein', {'missing': ... } )
        >>> k.add_residue( OBResidue.GetIdx(), 'protein', {'modified':  [133,...] } )
        >>> k.add_residue( OBResidue.GetIdx(), 'glyco', {'group':  [1,133,...] } )

        >>> k.get_chain('B')
        {0: {'protein', 'missing': ...}, }
'
        >>> k.types
        [ 'protein', 'ligand', 'dna', 'glyco', 'modifier', 'water' ]

        >>> k.tags
        [ 'missing','modified']

        # 'group'     : it should contain all members of the group (including self)

                        extracting any members of the group will return the whole
                        group, unless detach = True

        # 'modified'  : it should point to the modifying residue
        #

        it should become:
            0 : { 'type': 'protein'}

            1 : { 'type': 'glyco', 'info': {'mw. type whatever}, 'group':[1,2,3,4]}
            2 : { 'type': 'protein', 'modified':{1,} }
            3 : { 'type': 'glyco', 'group':[1,2,3,4]}
            4 : { 'type': 'glyco': 'group':[1,2,3,4]}

            5 : { 'type': 'protein', 'modified':{6,} }
            6 : { 'type': 'glyco', 'group':[5,6]}

    """

    def __init__(self, mol, data=None, debug=False):
        #DebugObj.__init__(self, debug)
        self.debug=debug
        self.mol = mol
        if data is None:
            data = {}
        self.data = data

    def add_residue(self, res_id, type_, info={}, replace=False):
        """ register a residue as a given type and associated info"""
        # print("========================================")
        # print("MESS CALLED ADD RESIDIUE", res_id, type_, replace)
        # pp(info)
        # print("========================================")
        if (res_id in self.data): # and (replace==False):
            # print("Collision! Residue[%d] already in the dataset!" % res_id)
            # raise KeyError
            data = self.data[res_id]
        else:
            res_obj = self.mol.GetResidue(res_id)
            name = res_obj.GetName()
            num = res_obj.GetNum()
            chain = res_obj.GetChain()

            data = {'type': [], 'name':name, 'num':num,'chain':chain, 'id':"%s:%s_%d" % (chain,name, num)}
        data.update(info)
        if not isinstance(type_,list):
            type_ = [type_]
        for t in type_:
            if not t in data['type']:
                # print("========>        APPENDING", t)
                data['type'].append(t)
        self.data[res_id] = data

    def get_residue(self, res_id):
        """ return data about a residue_id"""
        return self.data[res_id]

    def get_resname(self, res_id):
        """ return the residue name string CHAIN:RESnum"""
        resObj = self.mol.GetResidue(res_id)
        name = resObj.GetName()
        num = int(resObj.GetNumString())
        chain = resObj.GetChain()
        try:
            icode = resObj.GetInsertionCode().replace('\x00', '')
        except:
            icode = ""
        return "%s:%s%d%s" %(chain, name, num, icode)

    def extract_residue(self, res_id, complete=True):
        """ """
        res_list = [res_id]
        if not complete:
            return res_list
        d = self.data[res_id]
        if 'modified' in d:
            for modifier in d['modified']:
                res_list.append(modifier)
                mod_info = self.data[modifier]
                if 'group' in mod_info:
                    res_list += mod_info['group']
        return sorted(list(set(res_list)))

    def get_chain(self, chain_id):
        """ return all data associated with a chain"""
        data = {k:v for (k,v) in list(self.data.items()) if self.mol.GetResidue(k).GetChain() == chain_id}
        if len(data):
            return data
        return None

    def get_chains(self):
        """ """
        chains = []
        for r in list(self.data.keys()):
            chains.append(self.mol.GetResidue(r).GetChain())
        return list(set(chains))

    def change_data_type(self, type_=None, chain_id=None, res_id=None, unique=False):
        """ change the data type of a chain or a residue; by default the new type will be
            added to the existing ones, unless the `unique=True` flag is used
        """
        if type_ == None:
            raise TypeError("type_ argument must be specified")
        if chain_id is None and res_id is None:
            raise TypeError("Either chain_id or res_id must be specified")
        if chain_id is not None and res_id is not None:
            raise TypeError("Either chain_id or res_id must be specified")
        if not chain_id == None:
            for idx in self.data.keys():
                if self.data[idx]['chain'] == chain_id:
                    if unique:
                        self.data[idx]['type'] = [type_]
                    else:
                        self.data[idx]['type'].append(type_)
            return
        if not res_id == None:
            if unique:
                self.data[res_id]['type'] = [type_]
            else:
                self.data[res_id]['type'].append(type_)
            return

    def get_types(self):
        """ get a summary of all the types that have been registered """
        found = []
        for info in self.data.values():
            found.extend(info['type'])
        # print("FONDO", found)
        return list(set(found))

    def get_type(self, type_):
        """ retrieve data associated with a given type (protein, dna, ...)"""
        data = {res_id:info for (res_id, info) in list(self.data.items()) if type_ in info['type']}
        if len(data):
            return data
        else:
            return None

    def get_tags(self, full=False):
        """ retrieve the list of all residue tags contained in the data
            except for the required ones (name, num, chain), and type,
            unless full=True is used
        """
        tags = []
        for info in list(self.data.values()):
            tags += info
        tags = list(set(tags))
        #tags = list(set([info.keys() for info in self.data.values() ]))
        if full == True:
            return tags
        skip = ['info'] + [ 'type', 'group', 'target', 'target_id' ] + _required
        short = [ x for x in tags if not x in skip ]
        if len(short):
            return short
        return []

    def get_tag(self, tag):
        """ return data containing that property"""
        if tag in _required:
            print("Required tags cannot be used")
            raise Exception
        data = {res_id:info for (res_id, info) in list(self.data.items()) if tag in info }
        if len(data):
            return data
        else:
            return None

    def extract(self, types=None, single_res=None, exclude=None,
                exclude_seq_conflict=['EXPRESSION TAG', 'CLONING ARTIFACT'],
                detach=False, use_res_name=False, silent=True):
        """ used to extract information from the container
            if 'detach' is requested, then group information
            is ignored.

            # return protein and ligand types
            .extract(types=['protein','ligand'])

            # return residues 0,12,24
            .extract(singleResidue=[0,12,24])

            # extract all protein residues, except for 0,1,2
            .extract(types=['protein'], exclude=[0,1,2])

            # extract the whole dataset
            .extract(types='all')
        """
        if types == None:
            types = []
        elif types == 'all':
            types = self.get_types()
        if single_res == None:
            single_res = []
        if exclude == None:
            exclude = []
        # check that types are present
        if not silent:
            missing = []
            found_types = self.get_types()
            for t in types:
                if not t in found_types:
                    missing.append(t)
            if len(missing):
                msg = "Types not found in profile: %s" % ",".join(missing)
                raise KeyError(msg)
        accepted = []
        # extract types
        complete = not detach
        for t in types:
            type_data = self.get_type(t)
            if type_data == None:
                continue
            for k, v in list(type_data.items()):
                if not k in exclude:
                    keep = True
                    # check that the residue is not a sequence
                    # artifact like expression tag, cloning
                    # artifact, etc...
                    d = self.data[k]
                    if 'mutation' in d:
                        for e in exclude_seq_conflict:
                            if e in d['mutation']['description']:
                                keep = False
                                break
                    if keep:
                        accepted += self.extract_residue(k, complete)
        # extract single residues
        for s in single_res:
            accepted+= self.extract_residue(s, complete)
        accepted = list(set(accepted))
        if use_res_name == True:
            accepted = [ self.get_resname(x) for x in accepted ]
        return accepted



# profile.extract( types=['protein', 'dna', 'glyco'], single_res=[], exclude=[], detach=False)
#

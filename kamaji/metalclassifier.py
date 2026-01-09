class MetalClassifier:
    """ Prody-based class to classify coordinating metals
        and check if their geometry is compatible
        with catalytic activity
    """

    def __init__(self, kamajiprofile):
        """ """
        self.profile = kamajiprofile
        self.mol = self.profile.mol
        self.coordinatingMetals = chem.coordinatingMetals

    def canElement(self, string):
        """ generate the canonical element string, e.g.:
            ZN, zn -> Zn
        """
        buff = string[0].upper()
        if len(string) == 2:
            buff += string[1].lower()
        return buff

    def resInSelection(self, selection, unique=False):
        """ generate the list of residues in a selection; 
            optionally, only unique residues are returned
        """
        resList =  list(zip( selection.getChids(), selection.getResnames(), selection.getResnums() ))
        resList = [ ('%s:%s%s' % (x[0], x[1], x[2]) ) for x in resList ]
        if unique:
            return list(set(resList))
        return resList

    def processAll(self):
        """ """
        types = self.profile.get_types()
        if 'metal_catalytic' in types:
            salt = self.scanMetal(mode='ion')
            for s in salt:
                print("UPDATE PROFILE of SALT", s)
        if 'cofactor' in types:
            self.scanMetal(mode='cofactor')


        if not chain == None:
            chainList = [chain]
        else:
            chainList = list(self.interface.classes.keys())
        for chId in chainList:
            print("\n\nAnalyzing chain [%s] " % chId) 
            data = self.interface.classes[chId]
            if 'catalyticMetal' in list(data.keys()):
                salt = self.scanMetals(chId,data['catalyticMetal'], mode='ion')
                for kId in salt:
                    target = self.interface.classes[chId].setdefault('genericMetal', {})
                    target[kId] = self.interface.classes[chId]['catalyticMetal'].pop(kId)
            if 'cofactor' in list(data.keys()):
                self.scanMetals(chId,data['cofactor'], mode='cofactor')

    def get_res_info(self, selection):
        """ retrieve commonly used residue info from a selection"""
        res_names = selection.getResnames()
        res_nums = selection.getResnums()
        chIds = selection.getChids()
        elements = [ x.getElement() for x in selection ]
        elements = list(map(self.canElement, elements)) 
        indices = [ x.getIndex() for x in selection]
        res_id = []
        for i in range(len(indices)):
            res_id.append( "%s:%s%s" % (chIds[i], res_names[i], res_nums[i]))
        return {'num':res_nums, 'sel':selection, 'name':res_names, 'chains':chIds,
               'res_id':res_id, 'elements':elements, 'indices':indices}


    def scanMetalIon(self, res_id=None):
        """ examine metals identified as separate ions"""
        metals = self.coordinatingMetals
        if res_id == None:
            pass
        
        else:
            resList = [ res_id ]




        


    def scanMetalCofactor(self, res_id=None):
        """ examine metals that are part of a cofactor 
            [ XXX DO WE NEED THIS? XXX ]
        """
        pass

    def scanMetal(self, res_id=None, mode='ion'): # chId, dataSet, mode='ion'):
        """ examine metals either as separate ions (mode='ion')
            or as part of a cofactor, like heme (mode='cofactor')
        """
        assert mode in ['ion', 'cofactor']
        remove = []
        
        if mode == 'cofactor':
            cofactors = self.profile.get_type('cofactor')
            for c in cofactor:
                print("COFACT", c)
        return

            
            


        for kId, kData in list(dataSet.items()):
            if mode == 'cofactor':
                info = kData['info']
                if not 'metal' in info['setup']:
                    continue
                else:
                    print("\tFound cofactor with metal : [ %s : \"%s\" ]" % (info['id']['accuracy'],
                        info['name']))
            # get residue info
            res_num = kData['num']
            resSelection = self.mol.select('chain %s resnum %d' % (chId, res_num))
            res_info = self.get_res_info(resSelection)
            elementList = []
            indicesList = []
            # isolate known catalytic metals
            if mode == 'ion':
                elementList = [ res_info['elements'][0] ]
                indicesList = [ res_info['indices'][0] ]
            else:
                for i in range(len(res_info['indices'])):
                    e = res_info['elements'][i]
                    if e in metals:
                        elementList.append(e)
                        indicesList.append(res_info['indices'][i])
            for element, index in zip(elementList, indicesList):
                # get current metal info
                metalInfo = metals[element]
                dist = metalInfo['dist']
                coordContacts = metalInfo['contacts']
                coordElements = metalInfo['coordinatingElements']
                count, uniqueCount, contacts = self.getCoordinatingAtoms(mode=mode,
                                        chId=chId, res_num=res_num, index=index, 
                                        dist=dist, coordElements=coordElements)
                print("\tresidue [ %s ] " % res_info['res_id'][0])
                print("\telement [ %s ] : max.distance [ %2.3f ] ; accepted contacts [ %s ]" %  (element,
                        dist, ', '.join([ str(x) for x in list(coordContacts.keys())])))
                print("\tatom contacts   : %s" % (uniqueCount))
                print("\tatom non-unique : %s" % (count))
                print("\tresidues in contact: %d  [ %s ]" % (len(contacts),  ', '.join(list(contacts.keys()))))
                
                if uniqueCount in list(coordContacts.keys()):
                    # update metal info with coordination info
                    kData['coordinating'] = metals[element]['contacts'][uniqueCount]
                    kData['coordinating']['count'] = count
                    kData['coordinating']['uniqueCount'] = uniqueCount
                    kData['coordinating']['contacts'] = contacts
                    # debug
                    print("\t    *** found coordination geometry ***")
                    geom = coordContacts[uniqueCount]
                    for x,y in list(geom.items()):
                        print("\t    %s : %s" % (x,y))
                    print("\n")
                else:
                    print("\tPossible additive metal with %d contact(s)\n\n" % uniqueCount)
                    print("\t", list(contacts.keys()))                  
                    remove.append(kId)
        return remove


    def getCoordinatingAtoms(self, mode, chId=None, res_num=None, index=None, dist=None, coordElements=None):
        """ find atoms within coordination distance of the metal
            and count contact atoms; carboxy groups will be counted
            as one even in bi-dentate coordination
        """
        if mode == 'ion':
            selString = ('(exwithin %2.3f of chain %s resnum %d) '
                        'and (protein)' % (dist, chId, res_num) )
        elif mode == 'cofactor':    
            selString = ('(exwithin %2.3f of index %d) and '
                         '(protein or chain %s resnum %d )' % (dist, 
                         index, chId, res_num ) )
        else:
            print("getCoordinatingAtoms> ERROR MODE!")
            raise NameError
            
        shell = self.mol.select(selString)
        if shell == None:
            return (0, 0, {})
        res_info = self.get_res_info(shell)
        contacts = {}
        for i in range(len(res_info['indices'])):
            element = res_info['elements'][i]
            res_id = res_info['res_id'][i]
            index = res_info['indices'][i]
            if res_info['elements'][i] in coordElements:
                found = contacts.setdefault(res_id, {'indices':[], 
                    'elements' : []} )
                found['indices'].append(index)
                found['elements'].append(element)
        count, uniqueCount = self.countContacts(contacts)
        return (count, uniqueCount, contacts)

    def countContacts(self, contacts):
        """ count contacts checking that carboxy
            are considered only once
        """
        count = 0
        uniqueCount = 0
        for res_id, items in list(contacts.items()):
            firstCarboxy = True
            for i, element in enumerate(items['elements']):
                if not element == 'O':
                    count += 1
                    uniqueCount += 1
                    continue
                idx = items['indices'][i]
                if self.isCarboxy(idx):
                    if firstCarboxy:
                        count += 1
                        uniqueCount += 1
                        firstCarboxy = False
                    else:
                        count += 1
                else:
                    count += 1
                    uniqueCount += 1
        return count, uniqueCount
                        
                
    def isCarboxyCarbon(self, idx, exclude):
        """check that its a carbon in carboxyl group"""
        selString = 'bonded 1 to index %d'
        for a in self.mol.select(selString % idx):
            if not a.getIndex() == exclude:
                continue
            if not a.getElement() == 'O':
                continue
            return True
        return False


    def isCarboxy(self, index):
        """ check if the atom at index is an oxygen 
            in a carboxy atom
        """
        selString = 'bonded 1 to index %d'
        for a in self.mol.select(selString % index):
            if self.isCarboxyCarbon(a.getIndex(), index):
                return True
        return True


                    


    def averageCarboxy(self):
        pass




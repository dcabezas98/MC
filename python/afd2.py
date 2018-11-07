# coding=utf-8
import utils; from utils import rf,wf
import re

class Transition:
    """Models a transition in a dfa.



    """
    
    def __init__(self, sourceState, destState, label):
        """Initialize Transition object.

        Args:

            sourceState (str): the state from which this transition begins, e.g. 'q2'

            destState (str): the state in which this transition ends, e.g. 'q5'

            label (str): the scanned symbols for which this transition
                is followed, e.g. 'G'

          

        """
        self.sourceState = sourceState
        self.destState = destState
        self.label = label
       

    def __str__(self):
        return 'sourceState: %s destState: %s label: %s' %\
            (self.sourceState, self.destState, self.label)

    def __repr__(self):
        return self.__str__()
    
    # Static method to unify the labels of compatible transitions
    @staticmethod
    def unify(tList):
        """Unify the transitions in a list of compatible transitions.

        Transitions are compatible if they have the same source state,
        destination state. That is, they
        may differ in their label but nothing else. Sometimes it is
        convenient to unify compatible transitions by transforming
        them into a single transition with a new label that
        incorporates all of the input transitions. This method returns
        a single transition that unifies a given list of transitions,
        which must themselves all be compatible.

        Args:

            tList (list of Transition objects): the list of compatible
                transitions that will be unified.

        Returns:

            Transition: a single unified transition representing all
                transitions in tList.

        """
        assert len(tList)>0
        first = tList[0]
        unifiedTrans = Transition(first.sourceState, first.destState, None)
        for t in tList:
            assert unifiedTrans.isCompatible(t)
        labels = [t.label for t in tList]
        unifiedTrans.label = ''.join(labels)
        return unifiedTrans

    def isCompatible(self, other):
        """Determine whether this transition is compatible with the other transition.

        Transitions are compatible if they have the same source state,
        destination state, write symbol, and direction. That is, they
        may differ in their label but nothing else. This method
        returns True if this transition is compatible with the other
        transition, and False otherwise.

        Args:

            other (Transition): A Transition object to be compared
                with the calling Transition object.

        Returns:

            bool: True if this transition is compatible with the other
                transition, and False otherwise

        """
        
        return self.sourceState == other.sourceState \
            and self.destState == other.destState 
          

    def __eq__(self, other):
        if self is other: return True
        if other==None: return False
        if not isinstance(other, Transition): return False
        return self.sourceState == other.sourceState \
                and self.destState == other.destState \
                and self.label == other.label 

    def __ne__(self, other):
        return not self==other

    def __lt__ (self, other):
        return (self.sourceState, self.destState, self.label) \
            < \
            (other.sourceState, other.destState, other.label)

    def __gt__ (self, other):
        return other.__lt__(self)

    def getKeyForUnify(self):
        """Returns a key that can be used for collecting compatible transitions.

        See isCompatible() for a description of transition
        compatibility.

        Returns:

            (src, dest, writeSym, dir): a 4-tuple of str, consisting
                of the four attributes that determine with the
                transitions are compatible with each other.

        """
        return (self.sourceState, self.destState)

    @staticmethod
    def reassembleFromUnifyKey(key, label):
        """Re-create a transition from its compatibility key and label.

        See isCompatible() for a description of transition
        compatibility. See getKeyForUnify() for a description of
        compatibility keys. This method takes the 4-tuple
        compatibility key and returns a new Transition object
        compatible with that key and the given label.

        Args:

            key (4-tuple of str: (src, dest, writeSym, dir)): A
            compatibility key.

            label (str): the label for the new transition

        Returns:

            Transition: New transition object with the given label and
                compatibility key.

        """
        return Transition(key[0], key[1], label)






class afd:
    """Un autómata finito determinista.

    """
    
   
    # Un blanco se representa como guion bajo

    blank = '_'

    epsilon = '-'

    # Símbolos especiales 

    anySym = '~' # cualquier símbolo
    notSym = '!' # cualquier símbolo excepto el siguiente
    stateSeparator = '->'
    labelSeparator = ':'
    commentStart = '#'

    validSymbols = {c for c in
        r"""$'"%&()*+-./0123456789<>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}"""}
    
   
    
    def __init__(self, description = None, name = None, allowImplicitReject=False, verbose=True):
        """Initialize TuringMachine object.

        Args:

            description (str): Una cadena de caractereres con la descripcion del AFD

            tapeStr (str): La palabra de entrada

            name (str): Un nombre 

            allowImplicitReject (bool):
        """

        # Las transiciones son un diccionario con los estados como índice. Cada entrada contendrá 
	# la lista de transiciones que salen de ese estado 

        self.states = None
        self.final = None
        self.alphabet = None
        self.transitions = None

        

        # The name of the machine can be used for debugging and meaningful output.
        self.name = name
        self.allowImplicitReject = allowImplicitReject
        self.verbose = verbose

        
        if description:
            self.read(description)
        self.checkAllSymbolsValid()


    def splitTransition(self, line):
        """Given a line in a Turing machine description, split into transition components.

        Args:

            line (str): a line in a Turing machine description

        Returns:

            4-tuple of str (label, actions, sourceState, destState):
                where label, sourceState, destState are attributes as
                described in the documentation for the Transition
                class, and actions is a string containing the write
                symbol, if any, and the direction.

        """
        
        (source,rest) = line.split('->')
        (dest,label) = rest.split(':')
              
      
        return (source.strip(), dest.strip(), label.strip())

    def extractStates(self,line):

	
    	self.states = line.split()

    def extractFinal(self,line):

	
    	self.final = line.split()


    def extractTransition(self, line):
        """Given a line in a Turing machine description, return a new Transition object described by that   line.

        Args:

            line (str): a line in a Turing machine description

        Returns:

            Transition: a new Transition object

        """
        (source, dest, label) = \
                self.splitTransition(line)
       
        return Transition(source, dest, label)





    @staticmethod
    def stripComments(lines):
        """Strip comments from a Turing machine description.

        A comment is anything after a '#' chaacter on a given line.

        Args:

            lines (list of str): list of lines in the Turing machine
                description.

        Returns:

            list of str: The same list of Turing machine description
                lines with comments removed.

        """

        return [x.split(afd.commentStart)[0] for x in lines]

    
    def read(self, afdString):
        """Build the states and transitions of a Turing machine from an ASCII description.

        This method creates the self.transitions dictionary attribute,
        and populates it with the transitions in the given
        description. Nothing is returned. This method can also deal
        with building blocks, recursively reading descriptions of any
        building blocks encountered and adding them to the current
        machine.

        Args:

            tmString (str): Turing machine description.

        """

        self.transitions = dict()
        # split on newlines
        afdLines = afdString.split('\n')
        # strip comments
        afdLines = afd.stripComments(afdLines)

        
        
      

        self.extractStates(afdLines[0])
        self.extractFinal(afdLines[1])
        self.alphabet = afdLines[2].strip()

        afdLines = [x.strip() for x in afdLines[3:]]
 
     

        for line in afdLines:
            if len(line)>0:
                t = self.extractTransition(line)
                self.addTransition(t)


    def save(self,filename):
        descrip = self.write()
        wf(filename,descrip)

    def writeTransition(self, t):
        """Convert a transition into Turing machine description format.

        Args:

            t (Transition): the Transition object to be converted to
                description format

        Returns:

            str: description format of the transition t

        """

        components = [t.sourceState, afd.stateSeparator, t.destState,
                      afd.labelSeparator, ' ', t.label]
       
       
        return ''.join(components)
            
    def write(self):
        """Convert the current Turing machine into description format.

        Returns:

            str: description format of the current machine, suitable
                for storing in a .df file.

        """

        lines = []
        if self.transitions == None: 
            return '[No transitions]'
        line = ' '.join(self.states)
        lines.append(line)
        line = ' '.join(self.final)
        lines.append(line)
        lines.append(self.alphabet)
        for tList in self.transitions.values():
            for t in tList:
                line = self.writeTransition(t)
                lines.append(line)
        lines[3:].sort()
        return '\n'.join(lines)
    

    @staticmethod
    def labelMatchesSymbol(symbol, label):
        """Return True if the given symbol is valid for a transition with the given label.

        Usually, this will return True if symbol is one of the
        characters in label, but this method also handles certain
        special cases, such as the special symbol that matches any
        character, and the use of '!' for "not".

        Args:

            symbol (str): a single character

            label (str): the label attribute of a Transition. See the
                Transition documentation.

        Returns:

            bool: True if the symbol is valid for a transition
                with the given label, and False otherwise.

        """

        if afd.anySym == label:
            return True
        elif label[0] == afd.notSym:
            if symbol not in label[1:]: 
                return True
        elif symbol in label:
            return True
        else:
            return False
        

    # is t a valid transition?

    @staticmethod
    def isValidTransition(c,t):
        """Return True if t is a valid transition for the current configuration.

        Args:

            t (Transition): a Transition object

        Returns:

            bool: True if the current scanned symbol matches the label
                of t, meaning that t is a transition that can be
                followed from the current configuration.

        """

	
        
        return afd.labelMatchesSymbol(c, t.label)


    # get a list of the possible transitions from the given state
    def getTransitions(self, state):
        """Return a list of the possible transitions from the given state.

        This ignores the scanned symbol. The returned transitions are
        all the transitions that could ever be followed from the
        given state.

        Args:

            state (str): a state in the Turing machine

        Returns:

            list of Transition objects: A list containing all
                transitions whose source state is the given state
                parameter. This could be the empty list.

        """
        return self.transitions.get(state, [])


    def getValidTransitions(self,c,state):
        """Return a list of all valid transitions from the current configuration.

        This is a list of all transitions from the current state whose
        label matches the current scanned symbol.

        Returns:

            list of Transition objects: A list containing all valid
                transitions from the current configuration. This could
                be the empty list.

        """
        transitionList = self.getTransitions(state)
        ts = []
        for t in transitionList:
            if afd.isValidTransition(c,t):
                ts.append(t)
        return ts

    

    def run(self, input):
        """Apply the given transition to the current configuration.

        This implements one computational step of the Turing machine,
        following the given transition to its destination state,
        writing a symbol onto the tape if necessary, and moving the
        head if necessary.

        Args:

            t (Transition): the Transition to be followed. If t is
                None and implicit rejection is permitted, the machine
                will transition into the project state.

        """

        state = self.states[0]
        step = 0        

        for c in input:
           if self.verbose:
                print ("Step " + str(step) + " Current state " + state)
                print ("Reading " + c)
           trans = self.getValidTransitions(c,state)
           if len(trans) == 0:
                if self.allowImplicitReject:
                    return 'NO'
                    if self.verbose:
                        print ("Input Rejected, No defined transition")
                else:
                    return "E"
                    print ("Error, No defined transition")
           elif len(trans) ==1 :
                t = trans[0]
                state = t.destState
                if self.verbose:
                     print  ("Moving to state " + state)
                step = step + 1

           else:
                return "E"
                print  ("Error, Determinist Automata with more than a transition applicable")

        if state in self.final:
            if self.verbose:
                print  ("Input accepted. Final State")
            return 'SI'
        else:
            if self.verbose:
                print  ("Input rejected. Non Final State")
            return 'NO'


        



    def printTransitions(self):
        """Print the transitions of this machine"""
        for t in self.transitions.values():
            print(t)

    def addTransition(self, t):
        """Add a transition to this machine.

        Args:

            t (Transition): the Transition object to be added
        """
        if self.transitions == None:
            self.transitions = dict()
        transitionList = self.transitions.setdefault(t.sourceState, [])
        transitionList.append(t)

    def checkSymbolIsValid(self, t, c):
        """Check if a given symbol is permitted in Turing machines.

        Nothing is returned, but a WcbcException is raised if the
        symbol is invalid.

        Args:

            t (Transition): the Transition in which c is used.

            c (str): a single character which is the symbol to be
                checked

        """
        if c not in afd.validSymbols:
            message = '''***Error***: The symbol {0} (ASCII value {1}) is not permitted in Turing machine alphabets. The full transition containing this error is:\n{2}'''.format(c, ord(c), t)
            raise utils.WcbcException(message)
            

    def checkAllSymbolsValid(self):
        """Check that all symbols used in this machine are permitted.

        Nothing is returned, but a WcbcException is raised if a
        symbol is invalid.
        """
        if self.transitions:
            for tList in self.transitions.values():
                for t in tList:
                    label = t.label
                    if label==afd.anySym:
                        continue
                    elif label[0]==afd.notSym:
                        label = label[1:]
                    for c in label:
                        self.checkSymbolIsValid(t, c)


   
    def sortLabelChars(self, s):
        """Sort the characters in the transition label.

        We are given a line of a Turing machine description, and an
        equivalent line is returned. The only difference is that if
        there are multiple characters in the transition label, these
        characters are sorted into ascending order. This brings the
        description into a standard form that can be compared against
        other machines more easily.

        Args:

            s (str): a line of a Turing machine description

        Returns:

            str: a line equivalent to s but with transition label
                characters sorted

        """
        (prefix, label) = s.split(afd.labelSeparator, 1)
    
        sortedLabel = ''.join(sorted(label))
        return prefix + afd.labelSeparator + sortedLabel

    def standardizeDescription(self, d):
        """Return standardized version of the given Turing machine description.

        Args:

            d (str): A Turing machine description.

        Returns:

            str: a standardized version of d in which comments and
                whitespace and empty lines have been removed, the
                lines are sorted into lexicographical order, and
                labels have been sorted.

        """
        lines = d.split('\n')
        # remove comments
        lines = afd.stripComments(lines)
        # remove all whitespace (not just leading and trailing whitespace)
        lines = [re.sub(r'\s+', '', line) for line in lines]
        # remove empty lines
        lines = [x for x in lines if x!='']
        # sort characters within each label
        lines = [self.sortLabelChars(x) for x in lines]
        # sort
        lines.sort()
        # return single string with lines separated by newlines
        return '\n'.join(lines)
     


    def unifyTransitions(self):
        """Unify all transitions in this machine

        Transitions with the same source and destination state can be
        unified into a single transition with a longer label. This can
        be a useful way to simplify machine descriptions and to
        standardize them.

        """
        if self.transitions == None: return
        newTransitionsDict = dict()
        for state, tList in self.transitions.items():
            # key is tuple of all except label, value is list of labels
            tDict = dict()
            for t in tList:
                key = t.getKeyForUnify()
                labels = tDict.setdefault(key, [])
                labels.append(t.label)
            newTList = []
            for key, labels in tDict.items():
                labelStr = ''.join(labels)
                newTrans = Transition.reassembleFromUnifyKey(key, labelStr)
                newTList.append(newTrans)
            newTransitionsDict[state] = newTList
        self.transitions = newTransitionsDict

        
class RegularGrammar:

    epsilon='-'
    transition=' -> '

    def __init__(self, description=None):

        self.variables=None
        self.alphabet=None
        self.productions=None
        
        if(description):
            self.readFromFile(description)

            
    def readFromFile(self, inputfile):
        """
        Lee una gramática desde un fichero, con el formato especificado en el guión.
        """

        with open(inputfile) as f:

            self.variables = f.readline()[:-1].split(' ')
            self.alphabet = list(f.readline()[:-1])
            
            self.productions = dict()
            
            for v in self.variables:
                self.productions[v] = []
                
            """ Solo funciona para regulares

            for line in f:
                source = line[0]
                aux = line.split(RegularGrammar.transition)[-1]
                if aux[-1] in self.variables:
                    label, destiny = aux[0:-1], aux[-1]
                else:
                    label = aux
                    destiny = None
                
                self.productions[source].append((label,destiny))
            """
            
            for line in f:
                source, destiny = line.replace('\n','').split(RegularGrammar.transition)
                self.productions[source].append(destiny)
                
    
                
    def writeToFile(self, outputfile):
        """
        Escribe una gramática en un fichero, con el formato especificado en el guión.
        """

        f = open(outputfile,'w+')
        f.write(' '.join(self.variables)+'\n')
        f.write(''.join(self.alphabet)+'\n')

        for k in self.productions:
            for destiny in self.productions[k]:
                f.write(RegularGrammar.transition.join([k,destiny])+'\n')
        f.close()

    def checkLinearIzdaDcha(self):
        """
        Comprueba si la gramática es lineal por la izquierda o lineal por la derecha, devuelve un par de booleanos,
        el primero indica si es lineal por la izquierda y el segundo si lo es por la derecha.
        """
        
        linIzda = True
        linDcha = True
        
        for k in self.productions:
            
            if not linIzda and not linDcha:
                return False, False
                
            if k not in self.variables: # Ni siquiera es tipo 2
                return False, False

            for p in self.productions[k]:

                # Comprobamos si es lineal por la derecha
                if linDcha:
                    valido = False
                    if all((a in self.alphabet+[RegularGrammar.epsilon] for a in p)):
                        valido = True
                    else:
                        for var in self.variables:
                            if p[-len(var):] == var: # Debe tener una variable al final
                                if all((a in self.alphabet for a in p[:-len(var)])):
                                        valido = True
                    linDcha = valido

                # Comprobamos si es lineal por la izquierda
                if linIzda:
                    valido = False
                    if all((a in self.alphabet+[RegularGrammar.epsilon] for a in p)):
                        valido = True
                    else:
                        for var in self.variables:
                            if p[:len(var)] == var: # Debe tener una variable al final
                                if all((a in self.alphabet for a in p[len(var):])):
                                    valido = True

                    linIzda = valido
                    
        return linIzda, linDcha

    
    def checkRegular(self):
        """
        Comprueba si la gramática es regular, esto es, si es lineal por la izquierda o por la derecha
        """
        
        i, d = self.checkLinearIzdaDcha()
        return i or d

    
    
    def toAFND(self):
        """
        Genera un AFND a partir de una gramática lineal por la derecha usando el algoritmo visto en clase
        """

        i, d = self.checkLinearIzdaDcha()
        
        afnd = afd()

        afnd.states = self.variables.copy() + [RegularGrammar.epsilon]
        afnd.alphabet = ''.join(self.alphabet)
        afnd.final = [RegularGrammar.epsilon]
        
        if d:            
            for k in self.productions:
                for p in self.productions[k]:
                    afnd.addTransition(Transition(k,p,RegularGrammar.epsilon))
                    while p not in afnd.states and len(p)>1:
                        afnd.states.append(p)
                        source = p
                        label = p[0]
                        p = p[1:]
                        destiny = p
                        afnd.addTransition(Transition(source,destiny,label))
                    if len(p)==1:
                        if p not in afnd.states:                      
                            afnd.states.append(p)
                            afnd.addTransition(Transition(p,RegularGrammar.epsilon,p))
                            
        return afnd
                            
                
    def fromAFND(self, afnd):
        """
        Genera una gramática a partir de un autómata usando el algoritmo visto en clase
        """
        
        self.variables = afnd.states
        self.alphabet = afnd.alphabet

        self.productions = dict()
            
        for v in self.variables:
            self.productions[v] = []
                
        for q in afnd.states:
            l = afnd.transitions.get(q) # lista de transiciones desde q en automata

            for t in l:
                source, destiny = t.sourceState, t.label+t.destState
                self.productions[source].append(destiny)

        for q in afnd.final:
            self.productions[q].append(RegularGrammar.epsilon)

# see testCheckTM() in checkTuringMachine.py for more detailed tests
def testafd():
    for (filename, inString, solution) in [
            ('aut.af', '00111', 'SI'),
            ('aut.af', '00110', 'NO'),
            ]:
        aut = afd(rf(filename))
        val = aut.run(inString)
        utils.tprint('filename:', filename, 'inString:', inString, 'result:', val)
        assert val == solution

def testGrammar():
    
    print("Gramática lineal por la derecha")
    grd = RegularGrammar('grammar_dcha.gr')
    print(grd.checkLinearIzdaDcha())

    print("Gramática lineal por la izquierda")
    gri = RegularGrammar('grammar_izda.gr')
    print(gri.checkLinearIzdaDcha())
    
    print("Gramática lineal por la izquierda y por la derecha")
    grid = RegularGrammar('grammar_lin_id.gr')
    print(grid.checkLinearIzdaDcha())
    
    print("Gramática independiente del contexto")
    gr = RegularGrammar('grammar.gr')
    gr.writeToFile('grammar2.gr')
    print(gr.checkLinearIzdaDcha())

    print("Pasar de autómata a gramática")
    aut = afd(rf('aut.af'))
    gr_aut = RegularGrammar()
    gr_aut.fromAFND(aut)
    gr_aut.writeToFile('gr_aut.gr')
    print("Escrito en gr_aut.gr")

    print("Pasar de gramática (lineal por la derecha) a autómata")
    aut_gr = grd.toAFND()
    aut_gr.save('aut_gr.af')
    print("Escrito en aut_gr.af")

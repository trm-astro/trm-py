"""
Pure Python module to parse formulae

Needs much more work.
"""


class Formula(object):

    def __init__(self, form):
        """
        Initialises a formula from an input string.
        """
        pass

class Fnode(object):
    """
    Defines a node used to represent formulae.
    
    opvarnum  -- operator, variable or number as a string.
    ninput    -- number of inputs to this node.
    value     -- if ninput < 0, this node is the end of a branch and
                 has a value = value
    input     -- list of ninput Fnodes
    """

    def __init__(self, opvarnum, ninput, value):
        self.opvarnum = opvarnum
        self.ninput   = ninput
        self.value    = value
        self.inputs   = []

    def check(self):
        """
        Checks whether the Fnode is valid
        """
        if self.ninput:
            if self.ninput != len(self.inputs):
                raise Formula_Error('Number of inputs in Fnode does not match number expected')
            for inp in inputs:
                if not isinstance(inp, Fnode):
                    raise Formula_Error('At least one input to Fnode was not itself an Fnode')

    def add_input(self, fnode):
        """
        Adds an input to the Fnode
        """
        self.inputs.append(fnode)

class Formula_Error(exceptions.Exception):
    """For throwing exceptions from the Formula module"""
    def __init__(self, value):
        self.value = value
        
    def __str__(self):
        return repr(self.value)

    

"""Module containing classes and functions for the parsing of AMBER prmtop files"""

from itertools import izip_longest

format_strings = { 'a' : '%{0}s',
                   'E' : '%{0}e',
                   'I' : '%{0}d' }
data_type = { 'a' : str,
              'E' : float,
              'I' : int }

class MissingParamError ( Exception ):
    pass

def chunks ( seq, size, fillval = None ):
    "Returns an iterator that goes through seq in chunks of size"
    args = [ iter ( seq ) ] * size
    return izip_longest ( *args, fillvalue = fillval )

class section ():
    def __init__ ( self, body ):
        lines = body.split ( '\n' )
        self.title = lines[0].strip()
        pos = lines[1].find ( '(' )
        self.format = lines[1].strip()[ pos+1:-1 ]

        #prmtop files can have one of three formats denoted by a, E or I in the format line
        #these are tested and information from the format line stored
        for i in 'a', 'E', 'I':
            if i in self.format:
                self.char = i
        self.n, self.size = self.format.split( self.char )
        
        self.nchars = int ( self.size.split('.')[0] ) #account for strings with '.'
        self.format_string = format_strings[self.char].format ( self.size )

        #to get the values of the field reduce the lines back to a single string and remove line breaks
        #then iterate through in chunks of the correct size
        vals_txt = ''.join ( lines[2:] ).replace ( '\n', '' )
        self.vals = [ ''.join ( i ) for i in chunks ( vals_txt, self.nchars ) ]
        self.vals = map ( data_type[self.char], self.vals )
        
    def __str__ ( self ):
        val = '%FLAG {0:74s}'.format ( self.title )
        val += '\n%FORMAT{0:73s}\n'.format ( '(' + self.format + ')' )
        if self.vals == []:
            val += '\n'
            return val
        for i in chunks ( self.vals, int ( self.n ) ):
            for j in i:
                if j != None:
                    # val += '%{0}s'.format ( self.nchars ) % j
                    val += self.format_string % j
            val += '\n'
        # val += '\n'
        return val

class prmtop ():
    def __init__ ( self, filename ):
        with open ( filename ) as self.f:
            body = self.f.read()

        self.sections = map ( section, body.split ( "%FLAG" )[1:] )

        self.sec_dict = { i.title : i for i in self.sections[1:] }
        self.sec_dict[ 'VERSION' ] = body.split( "%FLAG" )[0]

    def __getitem__ ( self, val ):
        return self.sec_dict[val]

    def __iter__ ( self ):
        return iter ( self.sections )

    def __str__ ( self ):
        val = self[ 'VERSION' ]
        for sec in self.sections:
            val += str ( sec )
        return val

    def get_angle_by_at ( self, at1, at2, at3 ):
        for i in chunks ( self['ANGLES_INC_HYDROGEN'].vals + self['ANGLES_WITHOUT_HYDROGEN'].vals, 4 ):
            if sorted ( [at1, at2, at3] ) == sorted ( [ ptoi ( int ( j ) ) for j in i[:3] ] ):
                parm_set = int ( i[3] )
                return self['ANGLE_FORCE_CONSTANT'].vals[parm_set-1], self['ANGLE_EQUIL_VALUE'].vals[parm_set-1]
        else:
            raise MissingParamError ( "Unable to find angle associated with atoms %d %d %d" % ( at1, at2, at3 )  )

def ptoi ( p ):
    """Convert pointers used in the prmtop file to atom numbers"""
    return p / 3 + 1

def itop ( p ):
    """Convert atom numbers used in the prmtop file to pointers"""
    return ( p - 1 ) * 3


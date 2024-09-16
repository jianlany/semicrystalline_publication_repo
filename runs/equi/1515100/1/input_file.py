'''
  Class module for reading input files for the semicrystalline domain builder.
'''
class InputFile:
    ''' Reads arguments from the input file (same as input_file.h).
    Parses a human-readable input file for program options.  Each line in
    the input file defines an option.  Valid formats are:
        key
        key value
        # commented line
    '''
    def __init__(self, path):
        self.options = {'input_file': path}
        with open(path) as f:
            for line in f:

                line = line.split('#')[0].split()
                if len(line) == 0:
                    continue
                elif len(line) == 1:
                    self.options[line[0]] = True
                else:
                    self.options[line[0]] = ' '.join(line[1:]) 
    def __getattr__(self, key):
        ''' Allows this.key as a shortcut for this.options[key]. '''
        if key in self.options:
            return self.options[key]
        else:
            print 'Error: option', key, ' was not specified in input.'
            exit(1)

    def is_option_set(self, key):
        ''' Returns whether or not an option is set or defined. '''
        return key in self.options

if __name__ == '__main__':
    ''' This code is not for standalone mode, but you can test it by calling 
    it directly '''
    inputfile = InputFile('input.sc')
    for key in inputfile.options:
        print '{}: {}'.format(key, inputfile.options[key])

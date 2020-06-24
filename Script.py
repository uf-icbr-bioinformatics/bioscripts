### Utilities used by all bioscript programs

import sys
import os.path
import Utils

### Class to define subcommands in program

class Command():
    _cmd = ""
    __doc__ = "Use this to document commands."
    c = None                    # Database cursor, if needed

### Main Script class

class Script():
    name = ""
    version = "1.0"
    copyright = """(c) 2018, A. Riva, ICBR Bioinformatics Core, University of Florida
          L. Boatwright, ICBR Bioinformatics Core, University of Florida"""
    usagefun = None
    errorNames = {}
    errorMsg = {}
    errorCode = 1
    _commands = {}
    _commandNames = []
    _options = {}
    _arguments = []

    def __init__(self, name, version="1.0", usage=None, errors=[], copyright=None):
        """Errors should be a list of tuples: (code, name, message)."""
        self.name = name
        self.version = version
        self.usagefun = usage
        self.errorMsg = {}
        self.errorNames = {}
        self.defineErrors([('ERR', "Internal error", "Internal error: {}"),
                           ('NOFILE', "Missing file arguments", "Input file(s) not specified."),
                           ('BADFILE', "File not found", "Input file `{}' does not exist."),
                           ('BADINT', "Bad integer", "`{}' is not an integer number."),
                           ('BADFLOAT', "Bad float", "`{}' is not a floating point number.")])
        self.errorCode = 100
        self.defineErrors(errors)
        if copyright:
            self.copyright = copyright
        self.init()

    def init(self):
        pass

    def defineErrors(self, errors):
        for e in errors:
            self.errorNames[self.errorCode] = e[1]
            if len(e) == 3:
                self.errorMsg[self.errorCode] = e[2]
            else:
                self.errorMsg[self.errorCode] = e[1]
            setattr(self, e[0], self.errorCode)
            self.errorCode += 1
    
    def showErrors(self):
        for (code, name) in Utils.get_iterator(self.errorNames):
            sys.stderr.write("{}: {}\n".format(code, name))

    def usage(self, what=None):
        if what:
            self.usagefun(what)
        else:
            self.usagefun()
        sys.stderr.write(self.copyright + "\n")
        sys.exit(0)

    def getOptionValue(self, args, opts):
        """If one of the options in `opts' is present in `args', returns a tuple containing the actual argument
and the following one (if any). Otherwise, returns None."""
        nargs = len(args)
        for i in range(nargs):
            if args[i] in opts:
                if i < nargs-1:
                    return (args[i], args[i+1])
                else:
                    return (args[i],)
        return None

    def standardOpts(self, args):
        """Process the standard arguments. If any of them are found, this function does not return."""
        harg = self.getOptionValue(args, ['-h', '--help', '-E', '-v', '--version'])
        if harg:
            opt = harg[0]
            if opt == '-h' or opt == '--help':
                if len(harg) == 1:
                    self.usage()
                else:
                    self.usage(what=harg[1])
            elif opt == '-v' or opt == '--version':
                sys.stdout.write("{} version {}\n".format(self.name, self.version))
                sys.exit(0)
            elif opt == '-E':
                if len(harg) == 1:
                    self.showErrors()
                else:
                    errcode = int(harg[1])
                    if errcode in self.errorNames:
                        sys.stderr.write(self.errorNames[errcode] + "\n")
                        sys.exit(errcode)
                    else:
                        sys.stderr.write("Unknown error code {}\n".format(errcode))
                sys.exit(0)

    def errmsg(self, code, *args):
        if code in self.errorMsg:
            if len(args) > 0:
                msg = self.errorMsg[code].format(*args)
            else:
                msg = self.errorMsg[code]
            sys.stderr.write("Error: " + msg + "\n")
            sys.exit(code)
        else:
            sys.stderr.write("Unknown error code {}\n".format(code))
            sys.exit(0)
            
    # Standardized argument parsing

    def getOpt(self, opt):
        if opt in self._options:
            return self._options[opt]
        else:
            return None

    def getArgs(self):
        return self._arguments

    def parseArgs(self, args, specs):
        """Parse the supplied list of arguments `args' according to `specs', a comma-delimited
list of strings such as the following:

  `o'  - An option called -o taking no arguments, ie a switch (if present, its value will be True);
  `+o' - An option called -o taking a single argument;
  `+#o' - An option called -o taking a single argument, converted to integer;
  `*o' - An option called -o taking multiple arguments (stored as a list);
  `*#o' - An option called -o taking multiple arguments, converted to integers (stored as a list).

Option values are stored in the _options dictionary and can be retrieved with the getOpt() method.
All other command-line arguments are stored in the _arguments list and can be retrieved with the
getArgs() method.
"""
        wanted = self._parseArgspecs(specs)
        self._arguments = []
        self._options = {}
        for opt in wanted.keys():
            m = wanted[opt]
            if m == 0:
                self._options[opt] = False
            elif m == 1:
                self._options[opt] = None
            elif m == 2:
                self._options[opt] = []

        mode = None
        prev = ""
        for a in args:
            if a[0] == '-' and mode != 1:
                opt = a[1:]
                if opt in wanted:
                    mode = wanted[opt]
                    if mode == 0:
                        self._options[opt] = True
                        mode = None
                    else:
                        prev = opt
                else:
                    sys.stderr.write("Warning: unknown option `{}'.\n".format(opt))
            elif prev:
                if mode == 1:
                    self._options[prev] = a
                    mode = None
                    prev = ""
                elif mode == 2:
                    self._options[prev] = self.toInt(a)
                    mode = None
                    prev = ""
                elif mode == 3:
                    self._options[prev].append(a)
                elif mode == 4:
                    self._options[prev].append(self.toInt(a))
            else:
                self._arguments.append(a)

    def _parseArgspecs(self, specstring):
        specs = specstring.split(",")
        wanted = {}
        for s in specs:
            mode = 0
            if s[0] == "+":
                mode = 1
                s = s[1:]
                if s[0] == "#":
                    mode = 2
                    s = s[1:]
            elif s[0] =="*":
                mode = 3
                s = s[1:]
                if s[0] == "#":
                    mode = 4
                    s = s[1:]
            wanted[s] = mode
        return wanted

    def toInt(self, x, units=False):
        mult = 1
        if units:
            (x, mult) = Utils.decodeUnits(x)
        try:
            return int(x) * mult
        except ValueError:
            self.errmsg(self.BADINT, x)

    def toFloat(self, x):
        v = Utils.parseFraction(x)
        if v is None:
            try:
                return float(x)
            except ValueError:
                self.errmsg(self.BADFLOAT, x)
        else:
            return v
        
    def isFile(self, x, error=True):
        if x == "-":
            return x
        if os.path.isfile(x):
            return x
        elif error:
            self.errmsg(self.BADFILE, x)
        else:
            return False

    def addCommand(self, className):
        name = className._cmd
        self._commands[name] = className
        self._commandNames.append(name)
        self._commandNames.sort()

    def findCommand(self, cmd):
        if cmd in self._commands:
            return self._commands[cmd]
        else:
            return None

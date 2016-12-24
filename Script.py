### Utilities used by all bioscript programs

import sys
import os.path

class Script():
    name = ""
    version = "1.0"
    copyright = "(c) 2016, A. Riva, DiBiG, ICBR Bioinformatics, University of Florida"
    usagefun = None
    errorNames = {}
    errorMsg = {}
    errorCode = 1

    def __init__(self, name, version="1.0", usage=None, errors=[]):
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
        for (code, name) in self.errorNames.iteritems():
            sys.stderr.write("{}: {}\n".format(code, name))

    def usage(self):
        self.usagefun()
        sys.stderr.write(self.copyright + "\n")
        sys.exit(0)

    def standardOpts(self, args):
        if '-h' in args or '--help' in args:
            self.usage()
        if '-v' in args or '--version' in args:
            sys.stderr.write("{} version {}\n".format(self.name, self.version))
            sys.exit(0)
        if '-E' in args:
            p = args.index('-E')
            p += 1
            if p < len(args):
                errcode = int(args[p])
                if errcode in self.errorNames:
                    sys.stderr.write(self.errorNames[errcode] + "\n")
                    sys.exit(errcode)
                else:
                    sys.stderr.write("Unknown error code {}\n".format(errcode))
                    sys.exit(0)
            else:
                self.showErrors()
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

    def toInt(self, x):
        try:
            return int(x)
        except ValueError:
            self.errmsg(self.BADINT, x)

    def toFloat(self, x):
        try:
            return float(x)
        except ValueError:
            self.errmsg(self.BADFLOAT, x)
        
    def isFile(self, x):
        if os.path.isfile(x):
            return x
        else:
            self.errmsg(self.BADFILE, x)

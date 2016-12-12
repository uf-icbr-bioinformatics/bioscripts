### Utilities used by all bioscript programs

import sys

class Prog():
    name = ""
    version = "1.0"
    usagefun = None
    errorMsg = {}

    def __init__(self, name, version="1.0", usage=None, errors=[]):
        """Errors should be a list of tuples: (code, name, message)."""
        self.name = name
        self.version = version
        self.usagefun = usage
        self.errorMsg = {}
        for e in errors:
            self.errorMsg[e[0]] = e[2]
            setattr(self, e[1], e[0])

    def standardOpts(self, args):
        if '-h' in args or '--help' in args:
            self.usagefun()
            sys.exit(1)
        if '-v' in args or '--version' in args:
            sys.stderr.write("{} version {}\n".format(self.name, self.version))
            sys.exit(0)
        if '-E' in args:
            p = args.index('-E')
            p += 1
            if p < len(args):
                errcode = args[p]
                self.errmsg(int(errcode), True)

    def errmsg(self, code, exit=False):
        sys.stderr.write("Error: " + self.errorMsg[code] + "\n")
        if exit:
            sys.exit(code)
            

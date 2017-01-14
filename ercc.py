#!/usr/bin/env python

def parseLine(ln):
    ln = ln.rstrip("\r\n")
    pa = ln.split("\t")
    return pa

def posInList(l, value):
    try:
        return l.index(value)
    except ValueError:
        return None

def tableToFile(table, filename, hdr=None):
    with open(filename, "w") as out:
        if hdr:
            out.write("\t".join([str(x) for x in hdr]) + "\n")
        for row in table:
            out.write("\t".join([str(x) for x in row]) + "\n")

class ERCC():
    name = ""
    ERCCclass = ""
    concentration1 = 0
    concentration2 = 0
    fc = 0
    log2fc = 0

    def __init__(self, name="", ERCCclass="", concentration1=0, concentration2=0, fc=0, log2fc=0):
        self.name = name
        self.ERCCclass = ERCCclass
        self.concentration1 = concentration1
        self.concentration2 = concentration2
        self.fc = fc
        self.log2fc = log2fc

    def dump(self, stream):
        stream.write("ERCC(name='{}', ERCCclass='{}', concentration1={}, concentration2={}, fc={}, log2fc={})".format(
                self.name, self.ERCCclass, self.concentration1, self.concentration2, self.fc, self.log2fc))

    def init(self, fields):
        """Initialize an ERCC object with a parsed line from the ERCC-analysis file."""
        self.name = fields[1]
        self.ERCCclass = fields[2]
        self.concentration1 = float(fields[3])
        self.concentration2 = float(fields[4])
        self.fc = float(fields[5])
        self.log2fc = float(fields[6])

    def factor(self, mixA, mixB):
        """Returns the factor by which to multiply an ERCC using mix `mixA' to make it equivalent to one using mix `mixB'."""
        if mixA == mixB:
            return 1
        elif mixA == '1':
            return self.fc
        else:
            return 1/self.fc

class ERCCdb():
    entries = {
        'ERCC-00002': ERCC(name='ERCC-00002', ERCCclass='D', concentration1=15000.0, concentration2=30000.0, fc=0.5, log2fc=-1.0),
        'ERCC-00003': ERCC(name='ERCC-00003', ERCCclass='D', concentration1=937.5, concentration2=1875.0, fc=0.5, log2fc=-1.0),
        'ERCC-00004': ERCC(name='ERCC-00004', ERCCclass='A', concentration1=7500.0, concentration2=1875.0, fc=4.0, log2fc=2.0),
        'ERCC-00009': ERCC(name='ERCC-00009', ERCCclass='B', concentration1=937.5, concentration2=937.5, fc=1.0, log2fc=0.0),
        'ERCC-00012': ERCC(name='ERCC-00012', ERCCclass='C', concentration1=0.11444092, concentration2=0.17166138, fc=0.67, log2fc=-0.58),
        'ERCC-00013': ERCC(name='ERCC-00013', ERCCclass='D', concentration1=0.91552734, concentration2=1.83105469, fc=0.5, log2fc=-1.0),
        'ERCC-00014': ERCC(name='ERCC-00014', ERCCclass='D', concentration1=3.66210938, concentration2=7.32421875, fc=0.5, log2fc=-1.0),
        'ERCC-00016': ERCC(name='ERCC-00016', ERCCclass='C', concentration1=0.22888184, concentration2=0.34332275, fc=0.67, log2fc=-0.58),
        'ERCC-00017': ERCC(name='ERCC-00017', ERCCclass='A', concentration1=0.11444092, concentration2=0.02861023, fc=4.0, log2fc=2.0),
        'ERCC-00019': ERCC(name='ERCC-00019', ERCCclass='A', concentration1=29.296875, concentration2=7.32421875, fc=4.0, log2fc=2.0),
        'ERCC-00022': ERCC(name='ERCC-00022', ERCCclass='D', concentration1=234.375, concentration2=468.75, fc=0.5, log2fc=-1.0),
        'ERCC-00024': ERCC(name='ERCC-00024', ERCCclass='C', concentration1=0.22888184, concentration2=0.34332275, fc=0.67, log2fc=-0.58),
        'ERCC-00025': ERCC(name='ERCC-00025', ERCCclass='B', concentration1=58.59375, concentration2=58.59375, fc=1.0, log2fc=0.0),
        'ERCC-00028': ERCC(name='ERCC-00028', ERCCclass='A', concentration1=3.66210938, concentration2=0.91552734, fc=4.0, log2fc=2.0),
        'ERCC-00031': ERCC(name='ERCC-00031', ERCCclass='B', concentration1=1.83105469, concentration2=1.83105469, fc=1.0, log2fc=0.0),
        'ERCC-00033': ERCC(name='ERCC-00033', ERCCclass='A', concentration1=1.83105469, concentration2=0.45776367, fc=4.0, log2fc=2.0),
        'ERCC-00034': ERCC(name='ERCC-00034', ERCCclass='B', concentration1=7.32421875, concentration2=7.32421875, fc=1.0, log2fc=0.0),
        'ERCC-00035': ERCC(name='ERCC-00035', ERCCclass='B', concentration1=117.1875, concentration2=117.1875, fc=1.0, log2fc=0.0),
        'ERCC-00039': ERCC(name='ERCC-00039', ERCCclass='C', concentration1=3.66210938, concentration2=5.49316406, fc=0.67, log2fc=-0.58),
        'ERCC-00040': ERCC(name='ERCC-00040', ERCCclass='C', concentration1=0.91552734, concentration2=1.37329102, fc=0.67, log2fc=-0.58),
        'ERCC-00041': ERCC(name='ERCC-00041', ERCCclass='D', concentration1=0.22888184, concentration2=0.45776367, fc=0.5, log2fc=-1.0),
        'ERCC-00042': ERCC(name='ERCC-00042', ERCCclass='B', concentration1=468.75, concentration2=468.75, fc=1.0, log2fc=0.0),
        'ERCC-00043': ERCC(name='ERCC-00043', ERCCclass='D', concentration1=468.75, concentration2=937.5, fc=0.5, log2fc=-1.0),
        'ERCC-00044': ERCC(name='ERCC-00044', ERCCclass='C', concentration1=117.1875, concentration2=175.78125, fc=0.67, log2fc=-0.58),
        'ERCC-00046': ERCC(name='ERCC-00046', ERCCclass='D', concentration1=3750.0, concentration2=7500.0, fc=0.5, log2fc=-1.0),
        'ERCC-00048': ERCC(name='ERCC-00048', ERCCclass='D', concentration1=0.01430512, concentration2=0.02861023, fc=0.5, log2fc=-1.0),
        'ERCC-00051': ERCC(name='ERCC-00051', ERCCclass='B', concentration1=58.59375, concentration2=58.59375, fc=1.0, log2fc=0.0),
        'ERCC-00053': ERCC(name='ERCC-00053', ERCCclass='B', concentration1=29.296875, concentration2=29.296875, fc=1.0, log2fc=0.0),
        'ERCC-00054': ERCC(name='ERCC-00054', ERCCclass='C', concentration1=14.6484375, concentration2=21.9726563, fc=0.67, log2fc=-0.58),
        'ERCC-00057': ERCC(name='ERCC-00057', ERCCclass='C', concentration1=0.01430512, concentration2=0.02145767, fc=0.67, log2fc=-0.58),
        'ERCC-00058': ERCC(name='ERCC-00058', ERCCclass='C', concentration1=1.83105469, concentration2=2.74658203, fc=0.67, log2fc=-0.58),
        'ERCC-00059': ERCC(name='ERCC-00059', ERCCclass='D', concentration1=14.6484375, concentration2=29.296875, fc=0.5, log2fc=-1.0),
        'ERCC-00060': ERCC(name='ERCC-00060', ERCCclass='B', concentration1=234.375, concentration2=234.375, fc=1.0, log2fc=0.0),
        'ERCC-00061': ERCC(name='ERCC-00061', ERCCclass='D', concentration1=0.05722046, concentration2=0.11444092, fc=0.5, log2fc=-1.0),
        'ERCC-00062': ERCC(name='ERCC-00062', ERCCclass='A', concentration1=58.59375, concentration2=14.6484375, fc=4.0, log2fc=2.0),
        'ERCC-00067': ERCC(name='ERCC-00067', ERCCclass='B', concentration1=3.66210938, concentration2=3.66210938, fc=1.0, log2fc=0.0),
        'ERCC-00069': ERCC(name='ERCC-00069', ERCCclass='D', concentration1=1.83105469, concentration2=3.66210938, fc=0.5, log2fc=-1.0),
        'ERCC-00071': ERCC(name='ERCC-00071', ERCCclass='C', concentration1=58.59375, concentration2=87.890625, fc=0.67, log2fc=-0.58),
        'ERCC-00073': ERCC(name='ERCC-00073', ERCCclass='B', concentration1=0.91552734, concentration2=0.91552734, fc=1.0, log2fc=0.0),
        'ERCC-00074': ERCC(name='ERCC-00074', ERCCclass='C', concentration1=15000.0, concentration2=22500.0, fc=0.67, log2fc=-0.58),
        'ERCC-00075': ERCC(name='ERCC-00075', ERCCclass='B', concentration1=0.01430512, concentration2=0.01430512, fc=1.0, log2fc=0.0),
        'ERCC-00076': ERCC(name='ERCC-00076', ERCCclass='C', concentration1=234.375, concentration2=351.5625, fc=0.67, log2fc=-0.58),
        'ERCC-00077': ERCC(name='ERCC-00077', ERCCclass='D', concentration1=3.66210938, concentration2=7.32421875, fc=0.5, log2fc=-1.0),
        'ERCC-00078': ERCC(name='ERCC-00078', ERCCclass='D', concentration1=29.296875, concentration2=58.59375, fc=0.5, log2fc=-1.0),
        'ERCC-00079': ERCC(name='ERCC-00079', ERCCclass='D', concentration1=58.59375, concentration2=117.1875, fc=0.5, log2fc=-1.0),
        'ERCC-00081': ERCC(name='ERCC-00081', ERCCclass='D', concentration1=0.22888184, concentration2=0.45776367, fc=0.5, log2fc=-1.0),
        'ERCC-00083': ERCC(name='ERCC-00083', ERCCclass='A', concentration1=0.02861023, concentration2=0.00715256, fc=4.0, log2fc=2.0),
        'ERCC-00084': ERCC(name='ERCC-00084', ERCCclass='C', concentration1=29.296875, concentration2=43.9453125, fc=0.67, log2fc=-0.58),
        'ERCC-00085': ERCC(name='ERCC-00085', ERCCclass='A', concentration1=7.32421875, concentration2=1.83105469, fc=4.0, log2fc=2.0),
        'ERCC-00086': ERCC(name='ERCC-00086', ERCCclass='D', concentration1=0.11444092, concentration2=0.22888184, fc=0.5, log2fc=-1.0),
        'ERCC-00092': ERCC(name='ERCC-00092', ERCCclass='A', concentration1=234.375, concentration2=58.59375, fc=4.0, log2fc=2.0),
        'ERCC-00095': ERCC(name='ERCC-00095', ERCCclass='A', concentration1=117.1875, concentration2=29.296875, fc=4.0, log2fc=2.0),
        'ERCC-00096': ERCC(name='ERCC-00096', ERCCclass='B', concentration1=15000.0, concentration2=15000.0, fc=1.0, log2fc=0.0),
        'ERCC-00097': ERCC(name='ERCC-00097', ERCCclass='A', concentration1=0.45776367, concentration2=0.11444092, fc=4.0, log2fc=2.0),
        'ERCC-00098': ERCC(name='ERCC-00098', ERCCclass='C', concentration1=0.05722046, concentration2=0.08583069, fc=0.67, log2fc=-0.58),
        'ERCC-00099': ERCC(name='ERCC-00099', ERCCclass='C', concentration1=14.6484375, concentration2=21.9726563, fc=0.67, log2fc=-0.58),
        'ERCC-00104': ERCC(name='ERCC-00104', ERCCclass='B', concentration1=0.22888184, concentration2=0.22888184, fc=1.0, log2fc=0.0),
        'ERCC-00108': ERCC(name='ERCC-00108', ERCCclass='A', concentration1=937.5, concentration2=234.375, fc=4.0, log2fc=2.0),
        'ERCC-00109': ERCC(name='ERCC-00109', ERCCclass='B', concentration1=0.91552734, concentration2=0.91552734, fc=1.0, log2fc=0.0),
        'ERCC-00111': ERCC(name='ERCC-00111', ERCCclass='C', concentration1=468.75, concentration2=703.125, fc=0.67, log2fc=-0.58),
        'ERCC-00112': ERCC(name='ERCC-00112', ERCCclass='D', concentration1=117.1875, concentration2=234.375, fc=0.5, log2fc=-1.0),
        'ERCC-00113': ERCC(name='ERCC-00113', ERCCclass='C', concentration1=3750.0, concentration2=5625.0, fc=0.67, log2fc=-0.58),
        'ERCC-00116': ERCC(name='ERCC-00116', ERCCclass='A', concentration1=468.75, concentration2=117.1875, fc=4.0, log2fc=2.0),
        'ERCC-00117': ERCC(name='ERCC-00117', ERCCclass='B', concentration1=0.05722046, concentration2=0.05722046, fc=1.0, log2fc=0.0),
        'ERCC-00120': ERCC(name='ERCC-00120', ERCCclass='C', concentration1=0.91552734, concentration2=1.37329102, fc=0.67, log2fc=-0.58),
        'ERCC-00123': ERCC(name='ERCC-00123', ERCCclass='A', concentration1=0.22888184, concentration2=0.05722046, fc=4.0, log2fc=2.0),
        'ERCC-00126': ERCC(name='ERCC-00126', ERCCclass='B', concentration1=14.6484375, concentration2=14.6484375, fc=1.0, log2fc=0.0),
        'ERCC-00130': ERCC(name='ERCC-00130', ERCCclass='A', concentration1=30000.0, concentration2=7500.0, fc=4.0, log2fc=2.0),
        'ERCC-00131': ERCC(name='ERCC-00131', ERCCclass='A', concentration1=117.1875, concentration2=29.296875, fc=4.0, log2fc=2.0),
        'ERCC-00134': ERCC(name='ERCC-00134', ERCCclass='A', concentration1=1.83105469, concentration2=0.45776367, fc=4.0, log2fc=2.0),
        'ERCC-00136': ERCC(name='ERCC-00136', ERCCclass='A', concentration1=1875.0, concentration2=468.75, fc=4.0, log2fc=2.0),
        'ERCC-00137': ERCC(name='ERCC-00137', ERCCclass='D', concentration1=0.91552734, concentration2=1.83105469, fc=0.5, log2fc=-1.0),
        'ERCC-00138': ERCC(name='ERCC-00138', ERCCclass='B', concentration1=0.11444092, concentration2=0.11444092, fc=1.0, log2fc=0.0),
        'ERCC-00142': ERCC(name='ERCC-00142', ERCCclass='B', concentration1=0.22888184, concentration2=0.22888184, fc=1.0, log2fc=0.0),
        'ERCC-00143': ERCC(name='ERCC-00143', ERCCclass='C', concentration1=3.66210938, concentration2=5.49316406, fc=0.67, log2fc=-0.58),
        'ERCC-00144': ERCC(name='ERCC-00144', ERCCclass='A', concentration1=29.296875, concentration2=7.32421875, fc=4.0, log2fc=2.0),
        'ERCC-00145': ERCC(name='ERCC-00145', ERCCclass='C', concentration1=937.5, concentration2=1406.25, fc=0.67, log2fc=-0.58),
        'ERCC-00147': ERCC(name='ERCC-00147', ERCCclass='A', concentration1=0.91552734, concentration2=0.22888184, fc=4.0, log2fc=2.0),
        'ERCC-00148': ERCC(name='ERCC-00148', ERCCclass='B', concentration1=14.6484375, concentration2=14.6484375, fc=1.0, log2fc=0.0),
        'ERCC-00150': ERCC(name='ERCC-00150', ERCCclass='B', concentration1=3.66210938, concentration2=3.66210938, fc=1.0, log2fc=0.0),
        'ERCC-00154': ERCC(name='ERCC-00154', ERCCclass='A', concentration1=7.32421875, concentration2=1.83105469, fc=4.0, log2fc=2.0),
        'ERCC-00156': ERCC(name='ERCC-00156', ERCCclass='A', concentration1=0.45776367, concentration2=0.11444092, fc=4.0, log2fc=2.0),
        'ERCC-00157': ERCC(name='ERCC-00157', ERCCclass='C', concentration1=7.32421875, concentration2=10.9863281, fc=0.67, log2fc=-0.58),
        'ERCC-00158': ERCC(name='ERCC-00158', ERCCclass='B', concentration1=0.45776367, concentration2=0.45776367, fc=1.0, log2fc=0.0),
        'ERCC-00160': ERCC(name='ERCC-00160', ERCCclass='D', concentration1=7.32421875, concentration2=14.6484375, fc=0.5, log2fc=-1.0),
        'ERCC-00162': ERCC(name='ERCC-00162', ERCCclass='C', concentration1=58.59375, concentration2=87.890625, fc=0.67, log2fc=-0.58),
        'ERCC-00163': ERCC(name='ERCC-00163', ERCCclass='D', concentration1=14.6484375, concentration2=29.296875, fc=0.5, log2fc=-1.0),
        'ERCC-00164': ERCC(name='ERCC-00164', ERCCclass='C', concentration1=0.45776367, concentration2=0.68664551, fc=0.67, log2fc=-0.58),
        'ERCC-00165': ERCC(name='ERCC-00165', ERCCclass='D', concentration1=58.59375, concentration2=117.1875, fc=0.5, log2fc=-1.0),
        'ERCC-00168': ERCC(name='ERCC-00168', ERCCclass='D', concentration1=0.45776367, concentration2=0.91552734, fc=0.5, log2fc=-1.0),
        'ERCC-00170': ERCC(name='ERCC-00170', ERCCclass='A', concentration1=14.6484375, concentration2=3.66210938, fc=4.0, log2fc=2.0),
        'ERCC-00171': ERCC(name='ERCC-00171', ERCCclass='B', concentration1=3750.0, concentration2=3750.0, fc=1.0, log2fc=0.0)}

    def init(self, filename):
        """Initialize the ERCC database from an ERCC-analysis file."""
        self.entries = {}
        with open(filename, "r") as f:
            f.readline()
            for line in f:
                fields = parseLine(line)
                if len(fields) > 5:
                    e = ERCC()
                    e.init(fields)
                    self.entries[e.name] = e
                    
    def dump(self, stream):
        first = True
        stream.write("{")
        for k, e in self.entries.iteritems():
            if first:
                stream.write("'" + k + "': ")
                first = False
            else:
                stream.write(",\n'" + k + "': ")
            e.dump(stream)
        stream.write("}\n")

    def find(self, name):
        if name in self.entries:
            return self.entries[name]
        else:
            return None

    def makeERCCtable(self, filename, column, attribute='fc', filtercol=None, filtervalue=0, namecol=0):
        """Return a list of lists in which the first element is an ERCC name,
the second is the value of the `attribute' field for that ERCC, and the third
is the value in column `column' from `filename' for that ERCC. If `filtercol'
is specified, only consider rows from `filename' having that column greater than
`filtervalue'. Both `column' and `filtercol' can be specified either as numbers or
names."""
        result = []
        with open(filename, "r") as f:
            hdr = parseLine(f.readline())
            if type(column) == str:
                column = posInList(hdr, column) or column
            if type(filtercol) == str:
                filtercol = posInList(hdr, filtercol) or filtercol
            for line in f:
                fields = parseLine(line)
                ename = fields[namecol]
                ercc = self.find(ename)
                if ercc:
                    good = True
                    if filtercol:
                        if float(fields[filtercol]) < filtervalue:
                            good = False
                    if good:
                        result.append([ename, getattr(ercc, attribute), float(fields[column])])
        return result


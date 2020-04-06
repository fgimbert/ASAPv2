import os


class LocalFiles:
    """
    Set up local files for a project

    Parameters
    ----------
    - prj_info: dict
        project information (prj name, prefixes, ...)
    - project: str
        name of the project
    - homedir: str
        home directory
    - prjhome: str
        project home directory (homdir/project)
    - prefixes: list
        series of a prefix that sepcifies a contents of calculation
    - srcinp: dict{key: str}
        input files for QE on localhost
    - srcpp_inp: dict{key: str}
        input files for pp on localhost
    - srcpprism_inp: dict{key: str}
        input files for pprism on localhost
    - srcout: dict{key: str}
        output files for QE localhost
    - srcpp_out: dict{key: str}
        output files for pp on localhost
    - srcpprism_out: dict{key: str}
        output files for pprism on localhost
    - srcesm1: dict{key: str}
        ESM output files (chg, pot. profiles) on localhost
    - srcrism1: dict{key: str}
        RISM output files (chg, pot. density profiles) on localhost
    - src1drism: dict{key: str}
        1D-RISM outout files (RDF files for molecules) on localhost
    - srccube
        Gauusian cube file on localhost
    """
    def __init__(self, prj_info):
        self.project = prj_info['project']
        self.prefixes = prj_info['prefixes']
        self.homedir = prj_info['homedir']
        self.prjhome = self.homedir+'/'+self.project
        # initialize and set dict parameters
        self.srcinp = {}
        self.srcout = {}
        self.srccube = {}
        self.srcesm1 = {}
        self.srcrism1 = {}
        self.src1drism = {}
        self.srcpp_inp = {}
        self.srcpp_out = {}
        self.srcpprism_inp = {}
        self.srcpprism_out = {}
        for prefix in self.prefixes:
            self.srcinp[prefix] = self.prjhome+'/'+prefix+'.inp'
            self.srcout[prefix] = self.prjhome+'/'+prefix+'.out'
            self.srccube[prefix] = self.prjhome+'/'+prefix+'.cube'
            self.srcesm1[prefix] = self.prjhome+'/'+prefix+'.esm1'
            self.srcrism1[prefix] = self.prjhome+'/'+prefix+'.rism1'
            self.srcpp_inp[prefix] = self.prjhome+'/'+prefix+'_pp.inp'
            self.srcpp_out[prefix] = self.prjhome+'/'+prefix+'_pp.out'
            self.src1drism[prefix] = self.prjhome+'/'+prefix+'.1drism'
            self.srcpprism_inp[prefix] = self.prjhome+'/'+prefix+'_pprism.inp'
            self.srcpprism_out[prefix] = self.prjhome+'/'+prefix+'_pprism.out'
        # dummy parameters for QE input
        self.pseudo = 'PseudoDir'
        self.outdir = 'OutDir'
        self.maxseconds = 'MaxSeconds'
        # Check directory existence
        if os.path.exists(self.prjhome):
            print(f"Local working directory ({self.prjhome}) exists.")
            print(f"File list: {'  '.join(os.listdir(self.prjhome))}")
        else:
            print(f"{self.prjhome} does not exists")
            os.makedirs(self.prjhome)
            print(f"make {self.prjhome}")
        return

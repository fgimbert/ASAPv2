#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import paramiko
import sys
from scp import SCPClient
import datetime
from pyscript.remotefiles import RemoteFiles
from getpass import getpass


class RemoteJob(RemoteFiles):
    """
    Tools for remote job managing.

    Parameters
    ----------
    - prj_info: dict
        project information (prj name, prefixes, ...)
    - job_info: dict
        job information (number of MPI process, job duration, ...)
    - remote: str
        remote host name
    - queue: str
        job queue on remote host
    - current_job: str
        list of current running jobs
    - done_job: str
        list of finished jobs
    - srcsub: str
        job submitting script on localhost
    - srcshell: dict{key: str}
            job scripts on localhost
    - cputime: datetime valuable
        wall limit
    - node: int
        number of nodes
    - mpinum: int
        number of MPI processes
    - ompnum: int
        number of OpenMP threads
    - p_img: int
        number of images for QE
    - p_kpt: int
        number of pools for QE
    - p_tsk: int
        number of task parallelization
    - p_bnd: int
        number of band parallelization
    """
    def __init__(self, prj_info, job_info, remote, queue=None, savepass = None):
        super().__init__(prj_info, remote)
        
        if savepass is None:
            
            print("")
            print("Passphrase: (empty if no passphrase !)")
            self.passphrase = getpass()
            if self.passphrase == "":
                self.passphrase = None

            print("")
            print("Password: (empty if no password !)")
            self.password = getpass()
            if self.password == "":
                self.password = None
        else:
            self.passphrase = savepass[remote]["passphrase"]
            self.password = savepass[remote]["password"]
            
        self.kwargs = {}
        id_connect = {"username":self.user, "port":self.port, "passphrase":self.passphrase, "password":self.password, "key_filename":self.keyfile}
        for param in id_connect:
            if id_connect.get(param) != None:
                self.kwargs[param] = id_connect.get(param)
                
        # File path
        self.srcsub = self.prjhome+'/job_submit'
        self.currjob = self.prjhome+'/currentjob.txt'
        self.donejob = self.prjhome+'/donejob.txt'
        self.srcshell = {}
        for prefix in self.prefixes:
            self.srcshell[prefix] = self.prjhome+'/sub_'+prefix+'.sh'
        # Job information
        if queue is not None:
            self.queue = queue
        self.cputime = job_info['cputime']
        self.maxseconds = self.cputime.seconds-600
        self.node = job_info['node']
        self.mpinum = job_info['mpinum']
        self.totproc = job_info['node']*job_info['mpinum']
        self.ompnum = job_info['ompnum']
        self.jobid = {}
        self.p_img = job_info['para_img']
        self.p_kpt = job_info['para_kpt']
        self.p_tsk = job_info['para_tsk']
        self.p_bnd = job_info['para_bnd']
        # Replace to actual value
        list_script = []
        list_script.extend(arg_key for arg_key in self.script)
        for idx, line in enumerate(list_script):
            self.script[idx] = line.replace('QUEUE', self.queue) \
                                   .replace('_NODE_', str(self.node)) \
                                   .replace('TOTPROC', str(self.totproc)) \
                                   .replace('MPINUM', str(self.mpinum)) \
                                   .replace('OMPNUM', str(self.ompnum)) \
                                   .replace('CPUTIME', str(self.cputime)) \
                                   .replace('BINDIR', self.bindir) \
                                   .replace('KPT', str(self.p_kpt)) \
                                   .replace('BND', str(self.p_bnd)) \
                                   .replace('TSK', str(self.p_tsk))

    def connecthost(self):
        """
        Connect remote host and check existence of remote files
        """

        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())
              
            client.connect(self.host, **self.kwargs)
            stdin, stdout, stderr = client.exec_command("ls {0}".format(self.workdir))
            errs = stderr.read().decode('utf-8')
            if errs:
                print(f"Remote working directory ({self.workdir}) dose not exist.")
                stdin, stdout, stderr = client.exec_command("mkdir -p {0}".format(self.workdir))
                eerrs = stderr.read().decode('utf-8')
                if eerrs:
                    print(f"mkdir failed in connecthost: {self.workdir}")
                print(f"Make directory: {self.workdir}")
            else:
                print(f"Remote working directory ({self.workdir}) exists.")
                stdin, stdout, stderr = client.exec_command("ls {0}".format(self.workdir))
                eerrs = stderr.read().decode('utf-8')
                dir_client = stdout.read().decode('utf-8')
                if eerrs:
                    print(f"ls failed in connecthost: {self.workdir}")
                print(f"File list: {'  '.join(dir_client.split())}")
            client.close()
        return

    def copy_base(self, source_dir, target_dir):
        """
        Copy BASE data (wave function and potentials for restart) to current calculation prefix.
        """
        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())
            client.connect(self.host, username=self.user, port=self.port, passphrase=self.passphrase,
                           key_filename=self.keyfile)
            stdin, stdout, stderr = client.exec_command("ls {0}".format(target_dir))
            errs = stderr.read().decode('utf-8')
            if errs:
                print(f"Copy from BASE directory to: {target_dir}")
                stdin, stdout, stderr = client.exec_command("cp -r {0} {1}".format(source_dir, target_dir))
                eerrs = stderr.read().decode('utf-8')
                if eerrs:
                    print(f"COPY FAILED in copy_base: {eerrs}")
            else:
                print(f".save directory exists -> Remove and Copy BASE directory to: {target_dir}")
                stdin, stdout, stderr = client.exec_command("rm -fr {1}; cp -r {0} {1}".format(source_dir, target_dir))
                eerrs = stderr.read().decode('utf-8')
                if eerrs:
                    print(f"COPY FAILED in copy_base: {eerrs}")
            client.close()
        return

    def set_base(self, source_dir, target_dir):
        """
        Copy calculation data for restart to BASE directory.
        """
        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())
            client.connect(self.host, username=self.user, port=self.port, passphrase=self.passphrase,
                           key_filename=self.keyfile)
            stdin, stdout, stderr = client.exec_command("ls {0}".format(target_dir))
            errs = stderr.read().decode('utf-8')
            if errs:
                print(f"Set BASE directory to {target_dir}")
                stdin, stdout, stderr = client.exec_command("cp -r {0} {1}".format(source_dir, target_dir))
                eerrs = stderr.read().decode('utf-8')
                if eerrs:
                    print(f"COPY FAILED in set_base: {eerrs}")
            else:
                print(f"BASE directory exists -> Remove & Copy {target_dir}")
                stdin, stdout, stderr = client.exec_command("rm -fr {1}; cp -r {0} {1}".format(source_dir, target_dir))
                eerrs = stderr.read().decode('utf-8')
                if eerrs:
                    print(f"COPY FAILED in set_base: {eerrs}")
            client.close()
        return

    def put_file(self, srcfile):
        """
        Put local files (input files, job submit file and shell scripts) to remote working directory
        """
        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())
            client.connect(self.host, username=self.user, port=self.port, passphrase=self.passphrase,
                           key_filename=self.keyfile)
            with SCPClient(client.get_transport()) as scp:
                scp.put(srcfile, recursive=True, remote_path=self.workdir)
            client.close()
        return

    def get_file(self, tgtfile):
        """
        Get remote files (output files, esm1, rism1, 1drism files) to local project directory
        """
        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())
            client.connect(self.host, username=self.user, port=self.port, passphrase=self.passphrase,
                           key_filename=self.keyfile)
            with SCPClient(client.get_transport()) as scp:
                scp.get(tgtfile, recursive=True, local_path=self.prjhome)
            client.close()
        return

    def submit_job(self, subdirectory = ""):
        """
        Submit jobs and put JOBID into currentjob.txt
        """
        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())
            client.connect(self.host, **self.kwargs)
            stdin, stdout, stderr = client.exec_command("chmod +x {0}/job_submit".format(self.workdir + subdirectory))
            errs = stderr.read().decode('utf-8')
            if errs:
                print(f"chmod failed in submit_job: {self.workdir}")
            stdin, stdout, stderr = client.exec_command("cd {0}; ./job_submit".format(self.workdir + subdirectory))
            result = [line.split(self.jidsep)[self.jidcol] for line in stdout.readlines()]
            client.close()
        self.jobid = dict(zip(self.prefixes, result))
        for prefix in self.prefixes:
            print(f'Submit a job ({prefix}), Job ID = {self.jobid[prefix]}')

        now = datetime.datetime.now()
        try:
            with open(self.currjob, 'x', encoding='utf-8') as f:
                f.write('Host      Job ID    Project   Prefix    Status    Submit date   Starting data\n')
                for prefix in self.prefixes:
                    f.write(f'{self.name:<10}{self.jobid[prefix]:<10}{self.project:<10}{prefix:<10}queue     '
                            f'{now:%Y/%m/%d}\n')
        except FileExistsError:
            with open(self.currjob, 'a', encoding='utf-8') as f:
                for prefix in self.prefixes:
                    f.write(f'{self.name:<10}{self.jobid[prefix]:<10}{self.project:<10}{prefix:<10}queue     '
                            f'{now:%Y/%m/%d}\n')
        return

    def submit_files(self, prog, inpf, outf, basedir=None):
        """
        Generates job submit file (qsub XXX.sh) and shell scripts (mpirun -np XX ....)
        """
        f = open(self.srcsub, mode='w')
        for prefix in self.prefixes:
            f.write(self.qsub+" sub_"+prefix+".sh\n")
        f.close()

        for prefix in self.prefixes:
            f = open(self.srcshell[prefix], mode='w')
            list_script = []
            list_script.extend(arg_key for arg_key in self.script)
            if basedir is None:
                list_script.remove('COPY_BASE_DIR')
                for idx, line in enumerate(list_script):
                    newline = line.replace('PROG', prog) \
                                  .replace('PREFIX', prefix) \
                                  .replace('INPF', inpf[prefix]) \
                                  .replace('OUTF', outf[prefix])
                    f.write(newline)
            else:
                for idx, line in enumerate(list_script):
                    newline = line.replace('PROG', prog) \
                                  .replace('PREFIX', prefix) \
                                  .replace('INPF', inpf[prefix]) \
                                  .replace('OUTF', outf[prefix]) \
                                  .replace('COPY_BASE_DIR', f"cp -r {basedir} {self.outdir}/{prefix}.save\n")
                    f.write(newline)
            f.close()

        return

    def check_job(self, done):
        """
        Check the current running job on remote host. If a job is running the status is "1". If finished, it becomes 0.
        """
        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())
            client.connect(self.host, **self.kwargs)
            runid = []
            stdin, stdout, stderr = client.exec_command(self.qstat)
            result = stdout.readlines()
            client.close()
        for i, line in enumerate(result):
            if line == '\n':
                pass
            else:
                if line.split()[0] == self.jidchar:
                    for lline in result[i+self.jidskip:]:
                        runid.append(lline.split()[0])
                    break
        for prefix in self.prefixes:
            if self.jobid[prefix] not in runid:
                done[prefix] = True
        return

    def remote_qstat(self):
        """
        Excute 'qstat' on remote host.
        """
        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())
            client.connect(self.host, **self.kwargs)
            stdin, stdout, stderr = client.exec_command(self.qstat)
            result = stdout.readlines()
            for i in range(len(result)):
                print(result[i].rstrip('\n'))
            client.close()

    def check_queue(self):
        """
        Check the job status of a project.

        After submitting jobs, all JOBID is stored in the currentjob.txt.
        When 'check_queue' is called,
        * load the current running jobs from the currentjob.txt -> filejobs
        * do qstat and get current runnning jobs on a remote host -> runjobs
        * compare filejobs and runjobs, finished jobs go to donejob.txt, unfinished jobs go to currentjob.txt
        """
        # load running jobs from file and store it into dict.
        filejobs = {}
        try:
            with open(self.currjob, 'x', encoding='utf-8'):
                print(f"{self.currjob} does not exist.")
                sys.exit()
        except FileExistsError:
            with open(self.currjob, 'r', encoding='utf-8') as f:
                next(f)
                jobs = f.readlines()
                if jobs:
                    for job in jobs:
                        key = job.split()[0]+'_'+job.split()[1]  # host_JOBID
                        filejobs[key] = job.split()
                else:
                    print(f"{self.currjob} is empty. All jobs are done.")
                    sys.exit()
        # check current running jobs on remote host and store it into dict
        runjobs = {}
        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())
            client.connect(self.host, **self.kwargs)
            stdin, stdout, stderr = client.exec_command(self.qstat)
            result = stdout.readlines()
            for i, line in enumerate(result):
                if line == '\n':
                    pass
                else:
                    if line.split()[0] == self.jidchar:
                        for lline in result[i+self.jidskip:]:
                            key = lline.split()[0]  # JOBID
                            runjobs[key] = lline.split()[self.statcol]  # status only
                        break
        # Compare JOBID in runjobs with filejobs
        donejobs = {}
        for key in list(filejobs.keys()):
            if key.split('_')[0] == self.name:
                if key.split('_')[1] not in runjobs.keys():
                    donejobs[key] = filejobs[key]
                    del filejobs[key]
                else:
                    filejobs[key][4] = runjobs[key.split('_')[1]]

        now = datetime.datetime.now()
        # Store finished jobs in file
        try:
            with open(self.donejob, 'x', encoding='utf-8') as f:
                f.write('Host      Job ID    Project   Prefix    Submit date   End(checked) data\n')
                for key in donejobs.keys():
                    f.write(f'{donejobs[key][0]:<10}{donejobs[key][1]:<10}{donejobs[key][2]:<10}{donejobs[key][3]:<10}'
                            f'{donejobs[key][5]:<14}{now:%Y/%m/%d}\n')
        except FileExistsError:
            with open(self.donejob, 'a', encoding='utf-8') as f:
                for key in donejobs.keys():
                    f.write(f'{donejobs[key][0]:<10}{donejobs[key][1]:<10}{donejobs[key][2]:<10}{donejobs[key][3]:<10}'
                            f'{donejobs[key][5]:<14}{now:%Y/%m/%d}\n')

        # Store currently running jobs in file
        with open(self.currjob, 'w', encoding='utf-8') as f:
            f.write('Host      Job ID    Project   Prefix    Status    Submit date   Starting data\n')
            if filejobs.keys():
                print(f"Some running jobs are still runnning...")
                for key in filejobs.keys():
                    f.write(f'{filejobs[key][0]:<10}{filejobs[key][1]:<10}{filejobs[key][2]:<10}{filejobs[key][3]:<10}'
                            f'{filejobs[key][4]:<10}{filejobs[key][5]:<15}\n')
            else:
                print(f"All jobs are done!")
        return

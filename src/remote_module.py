from ipywidgets import widgets
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from IPython.display import Image, display, clear_output
# from pyscript.localfiles import LocalFiles
from traitlets import traitlets
# import SelectFilesWidget

# from SelectFilesWidget import *
from ipywidgets import Layout
from ipyfilechooser import FileChooser

import json
import paramiko
from scp import SCPClient

__author__ = "Florian Gimbert"
__copyright__ = "Copyright 2019, ASAP"
__version__ = "0.1"
__maintainer__ = "Florian Gimbert"
__email__ = "f-gimbert@nissan-arc.co.jp"
__status__ = "Development"
__date__ = "October 2019"


class LoadedButton(widgets.Button):
    """A button that can holds a value as a attribute."""

    def __init__(self, value=None, *args, **kwargs):
        super(LoadedButton, self).__init__(*args, **kwargs)
        # Create the value attribute.
        self.add_traits(value=traitlets.Any(value))


class Remote(object):

    def __init__(self, password=None, passphrase=None, gui=True):
        """Initialize remote class.
           Arguments
           ---------

        """

        self.clusters_list = []
        self.json_file = None

        self.passphrase = passphrase
        self.password = password
        self.host = None
        self.hostname = None
        self.port = None
        self.user = None
        self.keyfile = None
        self.workdir = None

        try:
            with open('remotehosts.json', mode='r', encoding='utf-8') as fd:
                self.json_file = json.load(fd)
                self.clusters_list = list(self.json_file.keys())
                # self.user = self.json_file[self.clusterslist_dropdown.value]['user']

        except FileNotFoundError:
            raise Exception("File 'remotehosts.json' does not exist")

        if gui:
            self.__create_remoteinfo_text()
            self.__add_cluster_button()
            self.__create_clusterslist_dropdown()
            self.__create_jobsdone_button()
            self.__create_jobsrunning_button()
            self.__create_boupdate_button()
            self.__create_bojob_text()
            self.__create_killjob_button()
            self.__create_restartjob_button()
            self.__create_outputjob_button()

            self.remote_button.on_click(self.addcluster_clicked)
            self.id_checkbox.observe(self.update_upload_button, 'value')
            self.clusterslist_dropdown.observe(self.update_informations, 'value')
            self.cluster_button.on_click(self.update_sshid)

        self.kwargs = {}
        self.id_connect = {"username": self.user, "port": self.port, "passphrase": self.passphrase,
                           "password": self.password, "key_filename": self.keyfile}

    def __create_remoteinfo_text(self):
        """Build the Text widget for Miller index  input"""
        self.remote_text = widgets.Text(
            placeholder='hostname',
            value=None,
            description='Hostname',
            disabled=False,
            style={'description_width': 'initial'}
        )

        self.host_text = widgets.Text(
            placeholder='host',
            value=None,
            description='Host',
            disabled=False,
            style={'description_width': 'initial'}
        )

        self.user_text = widgets.Text(
            placeholder='username',
            value=None,
            description='Username',
            disabled=False,
            style={'description_width': 'initial'}
        )

        self.port_text = widgets.Text(
            placeholder='port',
            value=None,
            description='Port',
            disabled=False,
            style={'description_width': 'initial'}
        )

        self.password_text = widgets.Password(description='Password:', value=None)
        self.passphrase_text = widgets.Password(description='Passphrase:', value=None)

        self.cluster_button = widgets.Button(description='Connect', disabled=False,
                                             tooltip='Connect to cluster information',
                                             icon='')

        self.cluster_validate = widgets.Valid(value=False, description='', readout='')


        self.id_checkbox = widgets.Checkbox(value=False,
                                            description='Select id key file',
                                            style={'description_width': 'initial'})

        # ***  upload id file button ***
        self.upload_button = FileChooser()
        self.upload_button.show_hidden = True


        #self.upload_button = widgets.FileUpload(
         #   description='id key file',
          #  disabled=False,
           # # accept='.txt',  # Accepted file extension e.g. '.txt', '.pdf', 'image/*', 'image/*,.pdf'
           # multiple=False,  # True to accept multiple files upload else False
           # layout=Layout(width='25%'),
           # button_style='success'
        #)

    def __read_remotehosts(self):

        # remote = RemoteFiles(prj_info=self.prj_info)
        print(self.clusters_list)

    def __add_cluster_button(self):
        """Build the button widget to add a remote cluster."""

        self.remote_button = LoadedButton(description='Add Remote Cluster', value=False, disabled=False)

    def __create_boupdate_button(self):
        """Build the button widget to print running jobs."""

        self.boupdate_button = widgets.Button(description='Update BO',
                                              disabled=False,
                                              button_style='success',
                                              tooltip='Update the list of BO running')

    def __create_jobsrunning_button(self):
        """Build the button widget to print running jobs."""

        self.jobsrunning_button = widgets.Button(description='Running jobs',
                                                 disabled=False,
                                                 button_style='success',
                                                 tooltip='Print list of running jobs')

    def __create_outputjob_button(self):
        """Build the button widget to print running jobs."""

        self.outputjob_button = widgets.Button(description='BO Ouput',
                                            disabled=False,
                                            button_style='success',
                                            tooltip='Show BO Output')

    def __create_restartjob_button(self):
        """Build the button widget to print running jobs."""

        self.restartjob_button = widgets.Button(description='Resart BO',
                                             disabled=False,
                                             button_style='warning',
                                             tooltip='Restart BO job')

    def __create_killjob_button(self):
        """Build the button widget to print running jobs."""

        self.killjob_button = widgets.Button(description='Kill BO',
                                             disabled=False,
                                             button_style='danger',
                                             tooltip='Kill BO job')

    def __create_bojob_text(self):
        """Build the button widget to print running jobs."""

        self.bojob_text = widgets.Text(description='BO id',
                                       disabled=False,
                                       placeholder='boid00')

    def __create_jobsdone_button(self):
        """Build the button widget to print finished jobs."""

        self.jobsdone_button = widgets.Button(description='Jobs done',
                                              disabled=False,
                                              button_style='warning',
                                              tooltip='Print list of finished jobs in the current session')

    def __create_clusterslist_dropdown(self):
        """Build the Dropdown widget for remote cluster selection"""

        self.clusterslist_dropdown = widgets.Dropdown(
            options=self.clusters_list,
            # placeholder='cluster',
            value=None,
            description='Available clusters',
            disabled=False,
            style={'description_width': 'initial'})

    def update_informations(self, *args):

        if self.clusterslist_dropdown is not None:
            try:
                with open('remotehosts.json', mode='r', encoding='utf-8') as fd:
                    self.json_file = json.load(fd)
                    # self.user = self.json_file[self.clusterslist_dropdown.value]['user']

            except FileNotFoundError:
                raise Exception("File 'remotehosts.json' does not exist")
            self.hostname = self.clusterslist_dropdown.value

            self.host = self.json_file[self.hostname]['host']
            self.port = self.json_file[self.hostname]['port']

            self.remote_text.value = self.hostname
            self.host_text.value = self.host

            if self.port == '':
                self.port = None

            self.user = self.json_file[self.hostname]['user']
            if self.port is None:
                self.port_text.value = ""
            else:
                self.port_text.value = self.port

            self.user_text.value = self.user
            if self.json_file[self.hostname]['keyfile'] is None:
                self.id_checkbox.value = False
            else:
                self.id_checkbox.value = True


    def update_upload_button(self, *args):

        if self.id_checkbox.value is True:
            self.upload_button.disabled = False
            self.upload_button.button_style = 'success'
        else:
            self.upload_button.disabled = True
            self.upload_button.button_style = 'warning'

    def update_sshid(self, *args):

        try:
            with open('remotehosts.json', mode='r', encoding='utf-8') as fd:
                self.json_file = json.load(fd)
                # self.user = self.json_file[self.clusterslist_dropdown.value]['user']

        except FileNotFoundError:
            raise Exception("File 'remotehosts.json' does not exist")

        if self.passphrase_text.value == '':
            self.passphrase = None
        else:
            self.passphrase = self.passphrase_text.value

        if self.password_text.value == '':
            self.password = None
        else:
            self.password = self.password_text.value

        

        self.hostname = self.clusterslist_dropdown.value

        self.host = self.json_file[self.hostname]['host']
        self.port = self.json_file[self.hostname]['port']
        if self.port == '':
            self.port = None

        self.user = self.json_file[self.hostname]['user']
        #self.workdir = self.json_file[self.hostname]['workdir']

        # self.keyfile = self.json_file[self.hostname]['keyfile']

        if self.id_checkbox.value is True:
            self.keyfile = self.upload_button.selected
        else:
            self.keyfile = None

        self.id_connect = {"username": self.user, "port": self.port, "passphrase": self.passphrase,
                           "password": self.password, "key_filename": self.keyfile}

        for param in self.id_connect:
            if self.id_connect.get(param) is not None:
                self.kwargs[param] = self.id_connect.get(param)


        cluster = {}
        kwargs2 = self.kwargs.copy()

        with open('../asapdata/cluster.json', mode='w', encoding='utf-8') as fp:

            cluster['host'] = self.host
            cluster['kwargs'] = kwargs2

            if 'password' in self.kwargs.keys():
                kwargs2['password'] = None
            if 'passphrase' in self.kwargs.keys():
                kwargs2['passphrase'] = None

            json.dump(cluster, fp)

    def read_remotehosts(self):

        # remote = RemoteFiles(prj_info=self.prj_info)
        print(self.clusters_list)

    def addcluster_clicked(self, b):

        if self.remote_text.value not in self.clusters_list:
            self.json_file[self.remote_text.value] = {'name': self.remote_text.value,
                                                      'host': self.host_text.value,
                                                      'user': self.user_text.value,
                                                      'port': self.port_text.value,
                                                      'keyfile': self.upload_button.selected}

            with open('remotehosts.json', mode='w', encoding='utf-8') as fp:
                json.dump(self.json_file, fp)

            self.clusters_list = list(self.json_file.keys())
            self.clusterslist_dropdown.options = self.clusters_list

    def command_ls(self, path=''):
        """
        Connect remote host and check remote files
        """
        # print(self.kwargs)
        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())

            client.connect(self.host, **self.kwargs)
            stdin, stdout, stderr = client.exec_command("ls {0}".format(path))
            errs = stderr.read().decode('utf-8')

            eerrs = stderr.read().decode('utf-8')
            dir_client = stdout.read().decode('utf-8')
            if eerrs:
                print(f"ls failed in connecthost:")
                return 1
            print(f"File list: {'  '.join(dir_client.split())}")
            client.close()

        return 0

    def command_qstat(self):
        """
        Connect remote host and check remote files
        """
        print(self.kwargs)
        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())

            client.connect(self.host, **self.kwargs)
            stdin, stdout, stderr = client.exec_command("qstat -u {0}".format(self.user))
            errs = stderr.read().decode('utf-8')

            eerrs = stderr.read().decode('utf-8')
            dir_client = stdout.read().decode('utf-8')
            if eerrs:
                print(f"ls failed in connecthost:")
            print(f"File list: {'  '.join(dir_client.split())}")
            client.close()

        return

    def put_file(self, path=None, files=None, host=None, kwargs=None):
        """
        Put local files (input files, job submit file and shell scripts) to remote working directory
        """

        if path is None:
            path = self.workdir
        if host is None:
            host = self.host
        if kwargs is None:
            kwargs = self.kwargs

        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())
            client.connect(host, timeout=30, **kwargs)

            with SCPClient(client.get_transport()) as scp:
                scp.put(files, recursive=True, remote_path=path)
            client.close()
        return

    def put_files(self, path=None, path_input=None, host=None, kwargs=None):
        """
        Put  local files from a repertory (input files, job submit file and shell scripts) to remote working directory
        """
        if path is None:
            path = self.workdir
        if host is None:
            host = self.host
        if kwargs is None:
            kwargs = self.kwargs

        import os
        arr = os.listdir(path_input)

        for file in arr:
            if file[0] is not '.':
                self.put_file(path=path, files=path_input + file, host=host, kwargs=kwargs)

        return


    def get_files(self, remote_path=None, local_path=None, host=None, kwargs=None):
        """
        Get  remote files from a repertory (input files, job submit file and shell scripts) to local working directory
        """
        if remote_path is None:
            path = self.workdir

        if host is None:
            host = self.host
        if kwargs is None:
            kwargs = self.kwargs

        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())
            client.connect(host, **kwargs)

            with SCPClient(client.get_transport()) as scp:
                scp.get(remote_path=remote_path, recursive=False, local_path=local_path)

        return

    def execute_command(self, command='ls', path=None, get_pty=False, host=None, kwargs=None):
        """
        Put local files (input files, job submit file and shell scripts) to remote working directory
        """
        if path is None:
            path = self.workdir
        if host is None:
            host = self.host
        if kwargs is None:
            kwargs = self.kwargs

        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())

            client.connect(host, **kwargs)
            stdin, stdout, stderr = client.exec_command("cd {0}/; {1}".format(path, command), get_pty=get_pty)
            #out = stdout.read().decode('utf-8')
            #print(out)
            #errs = stderr.read().decode('utf-8')
            #print(errs)


    def create_workdir(self, path=None, host=None, kwargs=None):
        """
        Put local files (input files, job submit file and shell scripts) to remote working directory
        """
        if path is None:
            path = self.workdir
        if host is None:
            host = self.host
        if kwargs is None:
            kwargs = self.kwargs

        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())

            client.connect(host, **kwargs)
            stdin, stdout, stderr = client.exec_command("ls {0}".format(path))
            errs = stderr.read().decode('utf-8')
            if errs:
                print(f"Remote working directory ({path}) does not exist.")
                stdin, stdout, stderr = client.exec_command("mkdir -p {0}".format(path))
                eerrs = stderr.read().decode('utf-8')
                if eerrs:
                    print(f"mkdir failed in connecthost: {path}")
                print(f"Make directory: {path}")
            else:
                print(f"Remote working directory ({path}) exists.")
                stdin, stdout, stderr = client.exec_command("ls {0}".format(path))
                eerrs = stderr.read().decode('utf-8')
                dir_client = stdout.read().decode('utf-8')
                if eerrs:
                    print(f"ls failed in connecthost: {path}")
                print(f"File list: {'  '.join(dir_client.split())}")

            client.close()
        return

    def execute_job(self, path=None, system=None, structure=None):
        """
        Execute job script inside workdir
        """
        print('Inside Execute Job')
        print('Path:', path)

        if path is None:
            path = self.workdir

        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy())

            client.connect(self.host, **self.kwargs)
            stdin, stdout, stderr = client.exec_command("ls {0}".format(path))
            errs = stderr.read().decode('utf-8')
            if errs:
                print(f"Remote working directory ({path}) does not exist.")
                stdin, stdout, stderr = client.exec_command("mkdir -p {0}".format(path))
                eerrs = stderr.read().decode('utf-8')
                if eerrs:
                    print(f"mkdir failed in connecthost: {path}")
                print(f"Make directory: {path}")
            else:
                print(f"Remote working directory ({path}) exists.")
                stdin, stdout, stderr = client.exec_command("chmod +x {0}/job_submit".format(path))
                stdin, stdout, stderr = client.exec_command("cd {0}/; ./job_submit".format(path))

                #stdin, stdout, stderr2 = client.exec_command("ls {0}".format(path))

                eerrs = stderr.read().decode('utf-8')
                dir_client = stdout.read().decode('utf-8')
                if eerrs:
                    print(f"ls failed in connecthost: {path}")
                    print(dir_client)
                else:
                    print(dir_client)
                    with open("running_jobs.txt", "a") as file:
                        file.write(system.rjust(10) + " | " + str(structure).rjust(10) + " | " + dir_client.split()[2].rjust(10) + '\n')

            client.close()
        return

    def update_jobs(self):

        import os

        if os.path.isfile('running_jobs.txt'):
            with open("running_jobs.txt", "r") as file:
                allines = file.readlines()

    def runningjobs_print(self):

        import os
        if os.path.isfile('running_jobs.txt'):
            print("System".rjust(10) + " | " + "Type".rjust(10) + " | " + "job-id".rjust(10))
            with open("running_jobs.txt", "r") as file:
                print(file.read())



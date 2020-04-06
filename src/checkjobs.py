import paramiko
import importlib
from src2 import optimizer_module, remote_module, viewer_mod, search_module, builder_module, vasp_module, data_module
import time
# Change in mymodule
importlib.reload(search_module)
importlib.reload(builder_module)
importlib.reload(optimizer_module)
importlib.reload(remote_module)
importlib.reload(vasp_module)
importlib.reload(data_module)
importlib.reload(viewer_mod)


def checkjobs(idjob, timecheck=10):
    # Execute VASP calculations inside test_i on enemat by ssh

    with paramiko.SSHClient() as client:
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.WarningPolicy)
        # Establish SSH connection
        client.connect(rjob.host, **rjob.kwargs)

        stdin, stdout, stderr = client.exec_command("qstat -r | grep -cw '{0}'".format(idjob))
        result = stdout.readlines()

    while int(result[0]) != 0:
        with paramiko.SSHClient() as client:
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.WarningPolicy)
            # Establish SSH connection
            client.connect(rjob.host, **rjob.kwargs)
            # stdin, stdout, stderr = client.exec_command("qstat -r | grep -cw {0}".format(prj_info['prefixes'][0]))
            stdin, stdout, stderr = client.exec_command("qstat -r | grep -cw '{0}'".format(idjob))
            result = stdout.readlines()

        time.sleep(timecheck*60)



    import os

    with open("jobs_finished.txt", "a") as file:
        file.write(system.rjust(10) + " | " + str(structure).rjust(10) + " | " + dir_client.split()[2].rjust(10) + '\n')


    if os.path.isfile('jobs_finished.txt'):
        print("System".rjust(10) + " | " + "Type".rjust(10) + " | " + "job-id".rjust(10))
        with open("jobs_finished.txt", "r") as file:
            print(file.read())
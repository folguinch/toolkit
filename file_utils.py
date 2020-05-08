import os
import subprocess

def write_txt(text, filename, mode='w'):
    with open(filename, mode) as out:
        out.write(text)

def rsync(origin, destination, *args):
    rsynccmd = 'rsync'
    flags = '-auvr'
    rsyncargs = [rsynccmd, flags] + list(args) + [origin, destination]
    subprocess.call(rsyncargs)
    
    if not os.path.exists(destination) and not os.path.isfile(destination):
        raise OSError('rsync failed')

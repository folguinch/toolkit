import smtplib, subprocess
from email.mime.text import MIMEText

from .decorators import timed

def sendMail(fromadd, toadd, subj, msg, smtp='phys.nthu.edu.tw', port=25):
    mess = MIMEText(msg)
    mess['Subject'] = subj
    mess['From'] = fromadd
    mess['To'] = toadd
    s = smtplib.SMTP(smtp, port)
    s.sendmail(fromadd, [toadd], mess.as_string())
    s.quit()

@timed
def run_mollie(dirname, logger=None, shell='csh', np=1, server=None,
        mpi='mpiexec', from_email=None, to_email=None):
    if shell=='csh':
        order = r"""csh -c 'cd %s; make; %s -n %i ./c.x'"""
    else:
        order = r"""cd %s; make; %s -n %i c.x"""
    if logger:
        logger.info('Running: %s', order % (dirname, mpi, np))
    if server:
        logger.info('Connecting to: %s', server)
        subp = subprocess.Popen(['ssh', server, order % (dirname,mpi,np)],
                                stderr=subprocess.PIPE, 
                                stdout=subprocess.PIPE)
    else:
        subp = subprocess.Popen(order % (dirname,mpi,np),
                                stderr=subprocess.PIPE, 
                                stdout=subprocess.PIPE,
                                shell=True)
    mollie_stdout, mollie_stderr = subp.communicate()
    if logger:
        logger.info('Mollie standard output:\n%s', mollie_stdout)
        logger.info('Mollie standard error:\n%s', mollie_stderr)

    # Send the stderr by email
    if mollie_stderr.startswith('Warning') and \
            len(mollie_stderr.splitlines())==0:
        pass
    elif from_email and to_email and mollie_stderr:
        msg = 'Order failed:\n%s\n\nMollie satandard error:\n%s' 
        msg = msg % (order % (dirname,mpi,np), mollie_stderr)
        subj = 'Mollie error'
        sendMail(from_email, to_email, subj, msg)

@timed
def run_casa(script, *args, **kwargs):
    casa = '/home/myso/opt/casa-release-4.7.2-el7/bin/casa'
    casa_opts = '--nologger --nogui --log2term'
    opts = ''
    for key, val in kwargs.items():
        if len(key)>1:
            opts += '--%s %s ' % (key, val)
        else:
            opts += '-%s %s' % (key, val)
    inp = ('%s '*len(args)).strip() % args
    cmd = '%s %s -c %s %s %s' % (casa, casa_opts, script, opts.strip(), inp)
    print cmd
    run_bash(cmd)

def run_bash(cmd, logger=None):
    subp = subprocess.Popen(cmd, stderr=subprocess.PIPE,
            stdout=subprocess.PIPE, shell=True)
    stdout, stderr = subp.communicate()
    if logger:
        logger.info(stderr)
    else:
        print stderr

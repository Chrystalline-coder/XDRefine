# -----------------------------------------------------------------------------
#  "THE BEER-WARE LICENSE" (Revision 43):
#  <lkrause@gwdg.de> wrote this file. As long as you retain this notice you
#  can do whatever you want with this stuff. If we meet some day, and you think
#  this is worth it, you can buy me a beer in return.
#  last change: 2016-09-05
#
#  README:
#  it's a python 'copy res2ins' routine
#  runs additional xd modules and calculates a 'multipole parameters to low-
#  angle data ratio'
#  there is a xdprop section to automatically set-up the .mas files, has to be
#  edited accordingly (default: CPSEARCH bond rmin  0.8 rmax  2.5)
#  
#  input: xd.hkl
#         xd01.inp
#         xd[01-99].mas
#  ----------------------------------------------------------------------------

low_cut = 0.5       # [  0.5] sinth/l value of data to be considered 'low order' in the calculation of low order to multipole parameter ratio (low2mp)
do_lsm = True       # [ True] do xd module every step
do_fft = False      # [False] do xd module every step
do_geo = False      # [False] do xd module every step
do_data2para = False # [ True] print d2p and low2mp after each step
do_final_fft = True # [ True] do xd module in last step
do_final_geo = True # [ True] do xd module in last step
do_final_fou = False # [ True] do xd module in last step
do_final_pro = False # [ True] do xd module in last step

import io,os,sys,re
import numpy as np
import shutil as sh
from glob import glob
from subprocess import Popen,PIPE,STDOUT
from datetime import datetime
curtime = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

# ################################################ #
#                     functions                    #
# ################################################ #
def where_to_start():
    master_files = sorted(glob('xd[0-9][0-9].mas'))
    complete_steps = sorted(glob('xd[0-9][0-9].res'))
    if len(complete_steps) == 0:
        return 1
    start_at = int(complete_steps[-1][2:4]) + 1
    if len(complete_steps) == len(master_files):
        start_at = 1
    while True:
        try:
            TEMP = raw_input('start from xd[{:>02}].mas? '.format(start_at)) or start_at
            if int(TEMP) <= 0 or int(TEMP) > len(master_files):
                raise ValueError
            return int(TEMP)
        except ValueError:
            continue

def run_prog(prog,XDName):
    p = Popen([prog,XDName],stdout=PIPE,stderr=STDOUT,bufsize=1)
    with p.stdout, open(logfile,'ab') as writer:
        for line in iter(p.stdout.readline, b''):
            sys.stdout.write(line)
            writer.write(line)
    p.wait()

def preprop():
    with open('xd.mas','r') as ofile:
        read_mas = ofile.readlines()
    replace = {'SELECT':'SELECT numdx *esd au verbose 1\n',
             #BzSe   'APPLY':'APPLY symm 2 translations 1 0 1 all\n', # 2nd half
             #BzSe calc cube at centre
             #BzSe'PROPERTY':'PROPERTY rho gradrho *d2rho nucpot core valence defden esp ef efg\n',
             #BzSe    'CUBE':'CUBE centre 2.14046 3.432679 7.015282 npts 90 stepsize 0.1\n',
             # critical point search
             'PROPERTY':'PROPERTY *rho gradrho d2rho nucpot core valence defden esp ef efg\n',
             'CPSEARCH':'CPSEARCH bond rmin  0.8 rmax  3.6\n'
              }
    active = False
    wfile = open('xd.mas','w')
    for line in read_mas:
        if 'MODULE XDPROP' in line:
            active = True
        if 'END XDPROP' in line:
            active = False
        if len(line.split()) > 0:
            if line.split()[0].strip('!') in replace and active:
                key = [key for key in replace if key in line and not replace[key] == None]
                wfile.write(replace.pop(key[0]))
            else:
                wfile.write(line)
    wfile.close()

def print_step_info(program,step):
    with io.open(logfile,'ab') as writer:
        sys.stdout.write('\n {:>^{w1}}\n >>>{:^{w2}}>>>\n >>>{:^{w2}}>>>\n >>>{:^{w2}}>>>\n'.format('>',' ','{}, step: {}'.format(program,step),' ',w1=66,w2=60))
        writer.write('\n {:>^{w1}}\n >>>{:^{w2}}>>>\n >>>{:^{w2}}>>>\n >>>{:^{w2}}>>>\n'.format('>',' ','{}, step: {}'.format(program,step),' ',w1=66,w2=60))

def data2para():
    with open('./xd_lsm.out') as ofile:
        read_lsmout = ofile.read().split('Residuals after final cycle')[0].split('Parameter / Variable Map')[1]
    exp_multipoles = re.compile('\s[MDQOH][01234][+-]*((?:\s+\d+)+)')# too lazy to figure out a regex that retrieves all \d+!
    find_mp_param = np.array(''.join(re.findall(exp_multipoles,read_lsmout)).split()).astype(int)# join+split for now
    find_mp_param[find_mp_param > 0] = 1
    mp_param = int(np.sum(find_mp_param))
    #r2 = re.search('R\{F\^2\}\s+=\s+(\d+\.\d+)\s+',read_lsmout).groups()[0]
    dp = round(float(re.search('Nref/Nv\s+=\s+(\d+\.\d+)\s+',read_lsmout).groups()[-1]),2)
    data = re.findall('\s+Included\s+in\s+the\s+refinement\s+=\s+(\d+)',read_lsmout)[0]
    para = int(re.search(' SCALE\s+(\d+)',read_lsmout).groups()[-1])
    q = re.search(' Rank of Q =\s+(\d+)',read_lsmout)
    if q:
        para = para - int(q.groups()[-1])
    with open('./xd.fco') as ofile:
        data_all_ar = np.genfromtxt(ofile,skip_header=26, usecols=(6,7))
    data_use_ar = data_all_ar[data_all_ar[:,1] == 0.]
    data_low_ar = data_use_ar[data_use_ar[:,0] <= low_cut]
    data_low = len(data_low_ar)
    if data_low > 0. and mp_param > 0.:
        low2mp = round(float(data_low)/mp_param,2)
    else:
        low2mp = 'None'
    return data,para,dp,data_low,mp_param,low2mp


# ################################################ #
#                      main                        #
# ################################################ #

with io.open('xd.hkl', 'rb') as xdhkl:
    XDName = xdhkl.readline().split()[0]
logfile = '{}_{}.log'.format(XDName,curtime)

if os.path.isfile('xd.inp') and not os.path.isfile('xd01.inp'):
    sh.copy('xd.inp', 'xd01.inp')

for i in range(where_to_start(),99):
    ii = '{:>02}'.format(i)

    #check end -> no more .mas files
    if os.path.isfile('xd{}.mas'.format(ii)):
        pass
    else:
        #last step done
        #recover latest results
        sh.copy('xd{:>02}.res'.format(i-1),'xd.res')
        #xdfft
        if do_final_fft:
            print_step_info('xdfft',i-1)
            run_prog('xdfft',XDName)
        #xdgeo
        if do_final_geo:
            print_step_info('xdgeom',i-1)
            run_prog('xdgeom',XDName)
        #xdfour
        if do_final_fou:
            print_step_info('xdfour',i-1)
            run_prog('xdfour',XDName)
        #xdprop
        if do_final_pro:
            print_step_info('xdprop',i-1)
            preprop()
            run_prog('xdprop',XDName)
        raise SystemExit

    #start
    sh.copy('xd{}.mas'.format(ii),'xd.mas')
    sh.copy('xd{}.inp'.format(ii),'xd.inp')

    #xdlsm
    if do_lsm:
        try:
            print_step_info('xdlsm',ii)
            run_prog('xdlsm',XDName)
            sh.copy('xd.res','xd{}.res'.format(ii))
            sh.copy('xd.inp','xd{}.inp'.format(ii))
            ##sh.move('xd.mas','xd{}.mas'.format(ii))
            ##sh.move('xd.mat','xd{}.mat'.format(ii))
            sh.copy('xd.cov','xd{}.cov'.format(ii))
            sh.copy('xd.fco','xd{}.fco'.format(ii))
            sh.copy('xd.fou','xd{}.fou'.format(ii))
            sh.copy('xd_lsm.out','xd{}_lsm.out'.format(ii))
            sh.copy('xd_lsm.cif','xd{}_lsm.cif'.format(ii))
            
            if do_data2para:
                data,para,dp,data_low,mp_param,low2mp = data2para()
                with io.open(logfile,'ab') as writer:
                    sys.stdout.write('\n {:^{w2}}\n'.format(' - d2p: {} ({}/{}) / low2mp: {} ({}/{}) -'.format(dp,data,para,low2mp,data_low,mp_param),w2=60))
                    writer.write('\n {:^{w2}}\n'.format(' - d2p: {} ({}/{}) / low2mp: {} ({}/{}) -'.format(dp,data,para,low2mp,data_low,mp_param),w2=60))

        except IOError:
            with open('xd_lsm.out') as file_lsm:
                lines = file_lsm.readlines()
            for line in reversed(lines):
                if '*' in line:
                    err_msg = line.strip() or "dumb"
                    break
                elif '-----' in line:
                    break
            with io.open(logfile,'ab') as writer:
                sys.stdout.write('\n {:>^{w1}}\n >>>{:^{w2}}>>>\n >>>{: ^{w2}}>>>\n {:>^{w1}}\n'.format('>','XD Error, step {}'.format(ii,XDName),err_msg,'>',w1=66,w2=60))
                writer.write('\n {:>^{w1}}\n >>>{:^{w2}}>>>\n >>>{: ^{w2}}>>>\n {:>^{w1}}\n'.format('>','XD Error, step {}'.format(ii,XDName),err_msg,'>',w1=66,w2=60))
            raise SystemExit

    #xdfft
    if do_fft:
        print_step_info('xdfft',ii)
        run_prog('xdfft',XDName)
        sh.move('xd_fft.out','xd{}_fft.out'.format(ii))
        sh.move('xd_fft.cif','xd{}_fft.cif'.format(ii))

    #xdgeom
    if do_geo:
        print_step_info('xdgeom',ii)
        run_prog('xdgeom',XDName)
        sh.move('xd_geo.cif','xd{}_geo.cif'.format(ii))
        sh.move('xd_geo.out','xd{}_geo.out'.format(ii))
        ##sh.copy('xd_geo.cry','xd{}_geo.cry'.format(ii))
        ##sh.copy('xd.tex','xd{}.tex'.format(ii))

    #remove xd.res, avoiding false positives
    sh.move('xd.res','xd{:>02}.inp'.format(i+1))

#!/usr/bin/python

import os
import sys
from subprocess import call, Popen, PIPE
from shutil import copy
import inspect
from socket import gethostname

main_dir = '/g/g17/dunn27/mgmol/'
cwd = os.getcwd()

this_file = inspect.stack()[0][1]

# Determine computer.
host = gethostname()
if 'vulcan' in host:
    queue = 'psmall'
    computer = 'bgq'
    walltime = '03:00:00'
    scratchdir = '/p/lscratchv/mgmolu/dunn27/mgmol'
    scratch = 'lscratchv'
    omp_num_threads = '4'
elif 'cab' in host:
    queue = 'pbatch'
    computer = 'pel'
    walltime = '01:00:00'
    scratchdir = '/p/lscratchd/dunn27/mgmol'
    scratch = 'lscratchd'
    omp_num_threads = '1'
else:
    raise SystemError('Unsupported computer.')

execfile(os.path.join(main_dir, 'util/submission/templates.py'))

prefix = cwd

if 'H2O_64' in prefix:
    example = 'H2O_64'
    md_template = md_template_H2O_64
    quench_template = quench_template_H2O_64
    params = {}
    params.update(H2O_64_params)

elif 'd144' in prefix or 'ShortSighted' in prefix or 'shortsighted' in prefix or 'equilibration' in prefix:
    example = 'ShortSighted'
    md_template = md_template_d144
    quench_template = quench_template_d144
    params = {}
    params.update(d144_params)

example_dir = os.path.join(main_dir, 'examples', example)

if 'd144' in prefix or 'ShortSighted' in prefix or 'shortsighted' in prefix or 'equilibration' in prefix:
    lrs = os.path.join(example_dir, 'lrs.in')
    copy(lrs, cwd)

coords = os.path.join(example_dir, 'coords.in')

copy(coords, cwd)

execfile(sys.argv[1])
num_runs = sys.argv[2]

def format_and_write(template, filename):

    line_iterator = iter(template.splitlines())
    
    with open(filename, 'wb') as f:

        category = ''

        for line in line_iterator:

            if line[0] == '[' and line[-1] == ']':
                try:
                    exec("%s = %s" % ('category', line.lstrip('[').rstrip(']')))
                except NameError:
                    pass

            elif '@' in line:
                attribute = line.split('=')[0]

                try:
                    value = category[attribute]
                except:
                    value = line.split('@')[1]
                    if len(value) == 0:
                        continue

                if type(value) == 'str':
                    value = value.lstrip("'").rstrip("'")
                else:
                    value = str(value)

                line = line.split('@')[0] + value

            f.write(line + '\n')

format_and_write(md_template, 'mgmol_md.cfg')
format_and_write(quench_template, 'mgmol_quench.cfg')

host = gethostname()

if 'cab' in host:
    params.update(cab_params)
elif 'vulcan' in host:
    params.update(vulcan_params)
else:
    raise RuntimeError('Unsupported machine.')

runfile_quench = runfile_quench_template.format(**params)
runfile_md = runfile_md_template.format(**params)

with open('runfile_quench', 'w') as f:
    f.write(runfile_quench)
with open('runfile_md', 'w') as f:
    f.write(runfile_md)


def jobs_log(job_dir):
    return os.path.join(job_dir, 'jobs.log')

def write_job_id(job_dir, job_id):
	with open(jobs_log(job_dir), 'a') as f:
		f.write(job_id + '\n')

def msub(job_dir, num_runs):
	print "Job Directory: "+job_dir
	
        log = jobs_log(job_dir)

	if os.path.exists(log):
		with open(log, 'r') as f:
			previous_job = f.readlines()[-1].lstrip().rstrip()
		output = Popen('msub -l depend=' + previous_job + ' ' + job_dir + '/runfile_md', stdout=PIPE, shell=True).communicate()[0]
	else:
		output = Popen('msub ' + job_dir + '/runfile_quench', stdout=PIPE, shell=True).communicate()[0]

	previous_job = output.lstrip().rstrip()
	write_job_id(job_dir, previous_job)

	for i in range(num_runs):
		print "Job #: "+str(i)
		output = Popen('msub -l depend=' + previous_job + ' ' + job_dir + '/runfile_md', stdout=PIPE, shell=True).communicate()[0]
		previous_job = output.lstrip().rstrip()
		write_job_id(job_dir, previous_job)

        run_string = 'msub -q pbatch -o mgmol.out -l nodes=1,walltime=00:05:00,depend=' + previous_job + ' -S $(which python) ' + this_file + ' ' + sys.argv[1] + ' ' + sys.argv[2]
        print run_string
	output = Popen(run_string, stdout=PIPE, shell=True).communicate()[0]
	previous_job = output.lstrip().rstrip()
	write_job_id(job_dir, previous_job)

msub(cwd, int(num_runs))

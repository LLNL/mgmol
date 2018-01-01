#!/usr/bin/env python

"""

prepare.py

Brief: Generic script for submitting jobs with a parameter sweep.

Usage: python prepare.py runs.param template.param

Authors: Guy Cohen, Yevgeny Bar Lev, Ian Dunn; Columbia University

"""

# Import modules.
from os import chdir, mkdir, getcwd,remove,path
from shutil import copy2
from glob import glob
import re
import itertools
import platform
from sys import argv
from subprocess import call
import textwrap
from datetime import timedelta,datetime
import inspect

# Data structure for parameters. 
class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

# Class for monitoring status of jobs.
class JobStatus:
    def __init__(self):
        self.total = 0
        self.done = 0
        self.queued = 0
        self.running = 0
        self.torun = 0
    def checkin(self, workdir):
        self.workdir = workdir
        if(not path.exists(cwd+'/'+workdir)):
            mkdir(cwd+'/'+workdir)
        self.total += 1
    def isDone(self):
        if path.exists(self.workdir+'/done'):
            self.done += 1
            return True 
    def isQueued(self):
        if path.exists(self.workdir+'/queued'):
            self.queued += 1
            return True
    def isRunning(self):
        if path.exists(self.workdir+'/running'):
            self.running += 1
            return True
    def setQueued(self):
        call(r'echo " " > '+cwd+'/'+self.workdir+'/queued', shell=True)
        self.torun += 1


def dir_name(prefix, pdict, secondary_keys=[], excluded_keys=[]):
    """A function for consistently naming data directories according to parameters."""
    def dict_to_name(d):
        
        items = d.items()

        for i in range(len(items)):

            items[i] = list(items[i])
            items[i][0] = items[i][0].split('.')[-1]
            items[i] = tuple(items[i])

        return '_'.join(map(lambda item : "%s%s" % item, items))

    main_items = pdict.copy()
    secondary_items = dict()
    for key in excluded_keys:
        del main_items[key]
    for key in secondary_keys:
        del main_items[key]
        secondary_items[key] = pdict[key]
    if len(secondary_keys)>0:
        return prefix+'_'+dict_to_name(main_items)+'/'+dict_to_name(secondary_items)
    else:
        return prefix+'_'+dict_to_name(main_items)


def param_string(pdict):
    """A function for creating a reduced parameter input file."""

    param_dict = {}

    for pair in pdict.items():

        genus = pair[0].split('.')[0]
        species = pair[0].split('.')[1]

        if genus not in param_dict.keys():

            param_dict[genus] = {}

        param_dict[genus][species] = pair[1]


    param_string = ""
   
    for pair in param_dict.items():

        param_string += pair[0]
        param_string += " = "
        param_string += str(pair[1])
        param_string += "\n"

    return param_string


if __name__ == "__main__":

    # Current working directory.
    cwd = getcwd()

    # Execute the parameter and template input files.
    if (len(argv)==1):
        this_file = inspect.stack()[0][1]
        this_dir = "/".join(this_file.split('/')[0:-1])
        execfile(path.join(this_dir, 'runs.param'))
        execfile(path.join(this_dir, 'template.param'))
    else:
        [execfile(argv[i]) for i in range(1,len(argv))]
    
    # Create an iterator which runs over all parameter sets:
    paramSets = []
    labels, terms = zip(*params.items())
    paramSets.append(itertools.product(*terms))
    paramIterator = itertools.chain.from_iterable(paramSets)
    
    # Label job array according to starting time.
    tag = datetime.now().strftime('%Y-%m-%d-%H%M%S')

    # Save list of jobs in job array in a file.
    job_list = open('array_job_list_'+str(tag),'w')
    js = JobStatus()

    # Sweep over all parameter sets, setup directories, and run jobs.
    for term in paramIterator:

        pdict = dict(zip(labels, term))    
        p = Struct(**pdict)#This makes it convenient to access the parameters as p.whatever...
        # name directories
        outdir = dir_name(prefix, pdict, secondary_keys=secondary_keys,excluded_keys=excluded_keys)
        # posting data    
        for valname,valfunc in post_data:
            pdict[valname] = valfunc(pdict)
    
        # check status 
        js.checkin(outdir)
        if js.isDone():
            print outdir, ' already done.'
        #elif(js.isQueued()):
            #print outdir, ' already queued.'
        elif(js.isRunning()):
            print outdir, ' already running.'
        else:
            js.setQueued()
            job_list.write(cwd + '/' + outdir + '\n')
    
            # formating template
            #p_contents = param_template.format(**pdict) #Commented out by Ian for flexible input files.
            p_contents = param_string(pdict)
            p_runfile = runfile_template.format(executable=executable)
            #if(not path.exists(cwd+'/'+outdir)):
                #mkdir(cwd+'/'+outdir)
            f = open(outdir+'/'+param_file_name, 'wb')
            f.write(p_contents)
            f.close()
            f = open(outdir+'/runfile', 'wb')
            f.write(p_runfile)
            f.close()
    
    job_list.close()
    run_array = run_array_template.format(num_jobs=js.torun,tag=tag)
    f = open(cwd+'/run_array_'+str(tag), 'wb')
    f.write(run_array)
    f.close()
    
    
    stat_content="""Job Statistics:
        Total:   {total}
        Done:    {done}
        Queued:  {queued}
        Running: {running}
    """
    print stat_content.format(total=js.total,done=js.done,queued=js.queued,running=js.running)
    if js.torun==0: print "no job to submit...exit"
    else:
        auth = raw_input("Submit {torun} jobs? (y/n) ".format(torun=js.torun))
        # Run the sbatch file
        if auth=='y':
        	call('bash run_array_'+str(tag),shell=True)
        #else:
            #remove('run_array_'+str(tag))
            #remove('array_job_list_'+str(tag))
    exit()

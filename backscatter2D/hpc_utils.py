import system_design
import numpy as np
import os, shutil
import re
import time

def setup_rotation_simulation_dirs(mcnp_input, num_theta, TR_cards = [3], axis='Z'):

    '''
    mcnp_input : name of mcnp input file to be copied and manipulated
                 for each rotation (type string)

    num_theta  : The number of steps to take during the scan rotation.
                If num_theta = 0 then there is only the one scan, at its
                original positon.
                If num_theta = 1 then there two total scans, one at the
                starting positon and one at 180 degrees from start.
                If num_theta = 3 then four total scans, one at every 90
                degrees

                etc....

                (type int)

    TR_cards   : list of the TR card names to be altered inside the mcnp
                 input file (type list, [type int] inside list)

    axis       : the axis which the rotation will occur.
                 Example: if axis Z then the X and Y axes will be altered
                          so that system is rotating about Z axis

                (type string)
                options:  'Z', 'X', 'Y'
    '''

    system_angles = np.linspace(0 ,
                                360.0*(1 - 1./(num_theta + 1)),
                                num_theta+1)

    pwd = os.getcwd()
    dir_list = []

    for degree in system_angles:

        path_name = os.path.join(pwd, str(int(degree))+'_degrees')

        mcnp_file = os.path.join(path_name, str(int(degree))+'_degrees.i')

        dir_list.append(path_name)

        try:
            os.mkdir(path_name)
            shutil.copyfile(mcnp_input, mcnp_file)
            os.chdir(path_name)

        except OSError:
            shutil.copyfile(mcnp_input, mcnp_file)
            os.chdir(path_name)

        f = open(mcnp_file, 'r')
        contents = f.read()
        f.close()

        for TR in TR_cards:

            regex = re.compile('\n\*TR' + str(TR) + '.*')

            original_TR = regex.findall(contents)[0]

            rotated_TR = __rotate_TR(original_TR, degree, axis)

            contents = contents.replace(original_TR, rotated_TR)
    

        contents = 'c\nc Imaging System Rotated by ' + str(degree) + \
                   ' degrees\nc\nc\n' + contents

        f = open(mcnp_file, 'w')
        f.write(contents)
        f.close()

        os.chdir(pwd)

    f = open('list_of_directories','w')

    for line in dir_list:
        f.write("%s\n" % line)

    f.close()

def __rotate_TR(TR_string, degree, axis):

    split_TR = TR_string.split(' ')

    name     = split_TR[0:4]
    
    x_axis   = np.array(map(float,   split_TR[4:7]))

    y_axis   = np.array(map(float,  split_TR[7:10]))

    z_axis   = np.array(map(float, split_TR[10:13]))

    end      = [split_TR[-1]]

    if axis == 'Z':

        new_x_axis = x_axis + np.array([degree, -degree, 0.])

        new_y_axis = y_axis + np.array([degree,  degree, 0.])

        return ' '.join(name + map(str, np.hstack([new_x_axis, 
                                                   new_y_axis, 
                                                   z_axis     ])) + end)

    if axis == 'X':

        new_y_axis = y_axis + np.array([0., degree, -degree])

        new_z_axis = z_axis + np.array([0., degree,  degree])

        return ' '.join(name + map(str, np.hstack([x_axis    , 
                                                   new_y_axis, 
                                                   new_z_axis ])) + end)


    if axis == 'Y':

        new_x_axis = x_axis + np.array([degree, 0., -degree])

        new_z_axis = z_axis + np.array([degree, 0.,  degree])

        return ' '.join(name + map(str, np.hstack([new_x_axis, 
                                                   y_axis    , 
                                                   new_z_axis ])) + end)

def add_slurm_files_to_dirs(job_sub, directory_list, num=None, resub=None, con_file='continue.i'):
    '''
    job_sub        = Name of .slurm file template to use to submit jobs with
                     on the HPC

    directory_list = Name of file that contains a listing of all the
                     directories that each simulation of interest is in

    num            = If you only want a certain index of directories from 
                     'directory_list' to add slurm files to, then specify 
                     the start and end index with num
                     Type: tuple (start int, end int)
 
    resub          = Type: string
                     options: 'resub' (If run didnt finish.
                                       Ex: because walltime expired)

                              'new_nps' (If you need to simulate
                                          more particles or want to stop at a
                                          new particile limit 'NPS')

    con_file       = If you choose to do 'new_nps' for resub, you need an
                     input file telling MCNP how many particles to simulate
    '''

    f = open(directory_list, 'r')
    dirs = [line.strip('\n') for line in f.readlines()]
    f.close()

    pwd = os.getcwd()

    if num == None:
        num = (0,len(dirs))

    if resub == None:
        for path in dirs[num[0]:num[1]]:

            job_name = path.split('/')[-1]

            os.system("cp %s %s" % (pwd+"/"+job_sub, path))

            os.system("sed -i 's/name.job/"+job_name+".job/g' "+\
                                                 path+"/"+job_sub)

            os.system("mv "+path+"/"+job_sub+" "+path+"/"+job_name+".slurm")

    elif resub == 'resub':
        for path in dirs[num[0]:num[1]]:

            job_name = path.split('/')[-1]

            os.system("cp %s %s" % (pwd+"/"+job_sub, path))

            os.system("sed -i 's/name.job/"+job_name+".job/g' "+\
                                                 path+"/"+job_sub)
            
            os.system("mv "+path+"/out1 "+path+"/old_out1")
            
            os.system("mv "+path+"/"+job_sub+" "+path+"/resub_"+\
                                                job_name+".slurm")
        
    elif resub == 'new_nps':
        for path in dirs[num[0]:num[1]]:
            
            job_name = path.split('/')[-1]
            
            os.system("cp %s %s" % (pwd+"/"+job_sub, path))
            
            os.system("cp %s %s" % (pwd+"/"+con_file, path))
            
            os.system("sed -i 's/name.job/"+job_name+".job/g' "+\
                                                 path+"/"+job_sub)
            
            os.system("mv "+path+"/out1 "+path+"/old_out1")
            
            os.system("mv "+path+"/"+job_sub+" "+path+"/new_nps_"+\
                                                 job_name+".slurm")
            
            os.system("mv "+path+"/"+con_file+" "+path+"/con"+job_name+".i")

def submit_slurm_jobs(directory_list, num=None, resub=None):
    '''
    directory_list = Name of file that contains a listing of all the
                     directories that each simulation of interest is in

    num            = If you only want a certain index of directories from 
                     'directory_list' to add slurm files to, then specify 
                     the start and end index with num
                     Type: tuple (start int, end int)

    resub          = Type: string
                     options: 'resub' (If run didnt finish. 
                                       Ex: because walltime expired)

                              'new_nps' (If you need to simulate
                                          more particles or want to stop at a
                                          new particile limit 'NPS')
    '''

    f = open(directory_list)
    dirs = [line.strip('\n') for line in f.readlines()]
    f.close()

    pwd = os.getcwd()

    if num == None:
        num = (0,len(dirs))

    if resub==None:
        for path in dirs[num[0]:num[1]]:

            job_name = path.split('/')[-1]

            os.chdir(path)

            #print "Submit Job "+job_name+".pbs from folder: "+os.getcwd()[-7:]
            os.system("sbatch --export=file_name="+job_name+".i "+job_name+".slurm")
            
            print('Submitting ' + path + '\n')
            
            time.sleep(1.0)

    if resub == 'resub':
        for path in dirs[num[0]:num[1]]:
            
            job_name = path.split('/')[-1]
            
            os.chdir(path)
            
            #print "Submit Job "+job_name+".pbs from folder: "+os.getcwd()[-7:]
            os.system("sbatch resub_"+job_name+".slurm")
           
            print('Resubmitting ' + path + '\n')

            time.sleep(1.0)
    
    if resub == 'new_nps':
        for path in dirs[num[0]:num[1]]:
            
            job_name = path.split('/')[-1]
            
            os.chdir(path)
            
            #print "Submit Job "+job_name+".pbs from folder: "+os.getcwd()[-7:]
            os.system("sbatch --export=file_name=con"+job_name+".i new_nps_"+job_name+".slurm")
            
            print('Continuing Simulation' + path + '\n')

            time.sleep(1.0)
        
    os.chdir(pwd)

def check_for_keyword_in_simulation_dirs(keyword=None, f_name=None, look_in='*.out'):
    '''
    keyword = The word you want to find in the 'look_in' file

    f_name  = Name of file you want to save the results in. If None,
              then it will just print the matches found inside the terminal

    look_in = the file (or type of file) you want to look inside of for the
              'keyword' inside all the simulation directories
    '''

    if f_name == None:
        os.system('grep "'+keyword+'" ./*/'+look_in)

    else:
        os.system('grep "'+keyword+'" -l ./*/'+look_in+' > temp')

        ft = open('temp','r')
        dirs = ft.read().split('\n')[:-1]
        ft.close()

        f = open(f_name,'w')
        pwd = os.getcwd()

        for path in dirs:
            job_name = path.split('/')[-2]
            f.write(pwd+'/'+job_name+'\n')

        f.close()
        os.system('rm temp')

def remove_files_in_simulation_dirs(directory_list, file_name, num=None):
    '''
    directory_list = Name of file that contains a listing of all the
                     directories that each simulation of interest is in

    file_name      = Name of file you wish to delete

    num            = If you only want a certain index of directories from 
                     'directory_list' to add slurm files to, then specify 
                     the start and end index with num
                     Type: tuple (start int, end int)
    '''
    f = open(directory_list)
    dirs = [line.strip('\n') for line in f.readlines()]
    f.close()

    pwd = os.getcwd()

    if num == None:
        num = (0,len(dirs))

    for path in dirs[num[0]:num[1]]:
        try:
            os.chdir(path)

            os.remove(file_name)

            os.chdir(pwd)

        except OSError:
            continue

def rename_files_in_simulation_dirs(directory_list, old_name, new_name, num=None):
    '''
    directory_list = Name of file that contains a listing of all the
                     directories that each simulation of interest is in

    old_name       = Current name of the file you wish to rename

    new_name       = New name you want to call the file

    num            = If you only want a certain index of directories from 
                     'directory_list' to add slurm files to, then specify 
                     the start and end index with num
                     Type: tuple (start int, end int)
    '''

    f = open(directory_list)
    dirs = [line.strip('\n') for line in f.readlines()]
    f.close()

    pwd = os.getcwd()

    if num == None:
        num = (0,len(dirs))

    for path in dirs[num[0]:num[1]]:
        try:
            os.chdir(path)

            os.rename(old_name, new_name)

            os.chdir(pwd)

        except OSError:
            continue

def print_last_dump_of_MCNP_simulations():

    pwd = os.getcwd()

    fdir = open('list_of_directories','r')

    dirs = fdir.readlines()

    fdir.close()

    dirs = [dir.strip('\n') for dir in dirs]

    for dir in dirs:

        os.chdir(dir)

        if os.path.isfile('out1'):

            f = open('out1', 'r')

            lines = f.readlines()

            f.close()

            print('\n'+ dir)

            print(lines[-3])

    os.chdir(pwd)

def collect_meshtallies_from_simulation_dirs(meshtally_file_name = 'meshtal',
                                        meshtally_file_dir = 'meshtally_files'):
    '''
    meshtally_file_name = Name of the meshtally files located in each simulation
                          direcory to be moved to the folder 
                          'meshtally_file_dir' and renamed

    meshtally_file_dir  = Name of directory to store all the collected meshtally
                          files
    '''

    f = open('list_of_directories', 'r')

    dirs = [line.strip('\n') for line in f.readlines()]

    f.close()

    pwd = os.getcwd()

    for path in dirs:

        mesh_name = path.split('/')[-1]

        file_location = os.path.join(path, meshtally_file_name)

        if os.path.isfile(file_location):

            new_file_location = os.path.join(pwd, meshtally_file_dir,
                                     '_'.join([mesh_name, 'meshtally']))

            try:
                shutil.copy2(file_location, new_file_location)

            except IOError:
                os.mkdir(meshtally_file_dir)

                shutil.copy2(file_location, new_file_location)

def recreate_list_of_directories_file(num_theta):
    '''
    num_theta  : The number of steps to take during the scan rotation.
                If num_theta = 0 then there is only the one scan, at its
                original positon.
                If num_theta = 1 then there two total scans, one at the
                starting positon and one at 180 degrees from start.
                If num_theta = 3 then four total scans, one at every 90
                degrees

                etc....

                (type int)
    '''

    system_angles = np.linspace(0 ,
                                360.0*(1 - 1./(num_theta + 1)),
                                num_theta+1)

    pwd = os.getcwd()

    dir_list = []

    for degree in system_angles:

        path_name = os.path.join(pwd, str(int(degree))+'_degrees')

        dir_list.append(path_name)
    
    f = open('list_of_directories','w')

    for line in dir_list:

        f.write("%s\n" % line)

    f.close()


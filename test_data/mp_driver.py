import multiprocessing as mp
import os
import subprocess
import shutil

global num_of_replications
global num_of_experiments
num_of_replications = 1
num_of_experiments = 100

def run_job(ijob):
    
    exp_num = ijob%num_of_experiments
    
    jobname = f"/lclscratch/murph/output_{ijob}"
    inputrecord = f"toplevel_outfiles/output_{ijob}.out"
    cwd = os.getcwd()
    base_dir = cwd
    suc_file = base_dir + '/compiled_data/output_' + str(exp_num) + '/succeeds.txt'
    fail_file = base_dir + '/compiled_data/output_' + str(exp_num) + '/fails.txt'
    pass_file = base_dir + '/compiled_data/output_' + str(exp_num) + '/passes.txt'
    
    # Initialize the output directory
    if not os.path.exists(base_dir + '/compiled_data/output_' + str(exp_num)):
        os.mkdir(base_dir + '/compiled_data/output_' + str(exp_num) )
        # Initialize the log files on how many attempts it takes to get the replications
        with open(fail_file, "w") as f:
            count = 0
            f.write(str(count) + "\n")
        with open(suc_file, "w") as f:
            count = 0
            f.write(str(count) + "\n")
        with open(pass_file, "w") as f:
            count = 0
            f.write(str(count) + "\n")

    # We wish to see if we have the correct number of replications.  If we do,
    # we do not need to run the algorithm again.
    #with open (suc_file, 'r') as f:
    #    num_sucs = int(f.read().strip())
        
    #if num_sucs >= num_of_replications:
    #    with open (pass_file, 'r') as f:
    #        chk_data = int(f.read().strip())
    #    with open(pass_file, "w") as f:
    #        count = chk_data + 1
    #        f.write(str(count) + "\n")
    
    # Otherwise, simulate dfnWorks to get an additional BT curve.
    # else:
    if True:
        with open(inputrecord, 'a') as f:
            param_names = "Starting job number:" + str(ijob) + "\n"
            f.write(param_names)
    
        ##########
        # The following try/except does not catch Segmentation Faults, which
        # appear to happen for MOST of our simulaitons.  We will use a series
        # of logfiles and analyses on the outputs to figure out when this happened.
        ##########
        try:
            cmd = f'python3.8 murph_dfn_v2.py {ijob} {num_of_experiments}  > outfiles/output_{ijob}.out'
            subprocess.call(cmd, shell=True)
            os.chdir(cwd)
            completed = 'true'
            print(f'makde it to the end of the job statement for job id: {ijob}')
        except:
            print(f'error with job id : {ijob}')
            os.chdir(cwd)
            completed = 'false'
    
        #completed = 'true'
        with open(inputrecord, 'a') as f:
            param_names = "Job completed:" + completed + "\n"
            f.write(param_names)    
            param_names = "Exiting job number:" + str(ijob)
            f.write(param_names)
    
        # From here, gather the important data values from the output, 
        # then delete the rest.
        print(f'transfering relevant output files for job {ijob} to a single location...')
        location = jobname + '/dfnTrans_output/partime.dat'
        destination = base_dir + '/compiled_data/output_' + str(exp_num) + '/'
        
        # We want to have at least 10% breakthrough.  If we don't have that, we
        # will consider this simulation a failure.
        bt_flag = False
        if os.path.exists(location):
            with open(location, 'r') as file:
                li = file.readlines()
                total_lines = len(li)
            # Check to see if we got at least 90% BT:
            if total_lines > 90002:
                shutil.move(location, destination)
                new_name = destination + 'partime' + str(ijob) + '.dat'
                old_name = destination + 'partime.dat'
                os.rename(old_name,new_name)
                
                with open (suc_file, 'r') as f:
                    chk_data = int(f.read().strip())
                with open(suc_file, "w") as f:
                    count = chk_data + 1
                    f.write(str(count) + "\n")
                flag = False
            else:
                bt_flag = True
                flag = True
        else:
            flag = True
    
        # In the case of a failure, increment the failure log.
        if flag:
            with open (fail_file, 'r') as f:
                chk_data = int(f.read().strip())
            with open(fail_file, "w") as f:
                count = chk_data + 1
                f.write(str(count) + "\n")
            with open(inputrecord, 'a') as f:
                param_names = "If (" + str(bt_flag) + ") then we didn't get 90k BT, otherwise the file was not written for some reason."
                f.write(param_names)
        # Otherwise, move over the additional files that we want.
        else:
            location = jobname + '/output_' + str(ijob) + '_effective_perm.txt'
            destination = base_dir + '/compiled_data/output_' + str(exp_num) +'/'
            if os.path.exists(location):
                shutil.move(location, destination)
                with open(inputrecord, 'a') as f:
                    param_names = "Moved" + location + "to " + destination
                    f.write(param_names)
    
            location = jobname + '/SA_params.txt'
            destination = base_dir + '/compiled_data/output_' + str(exp_num) + '/'
            if os.path.exists(location):
                shutil.move(location, destination)
                new_name = destination + 'SA_params_' + str(ijob) + '.dat'
                old_name = destination + 'SA_params.txt'
                os.rename(old_name,new_name)
                with open(inputrecord, 'a') as f:
                    param_names = "Moved" + location + "to " + destination + " and renamed it" + new_name
                    f.write(param_names)
                
            location = jobname + '/partime.dat'
            destination = base_dir + '/compiled_data/output_' + str(exp_num) + '/'
            if os.path.exists(location):
                shutil.move(location, destination)
                new_name = destination + 'graph_partime_' + str(ijob) + '.dat'
                old_name = destination + 'partime.dat'
                os.rename(old_name,new_name)
                with open(inputrecord, 'a') as f:
                    param_names = "Moved" + location + "to " + destination + " and renamed it" + new_name
                    f.write(param_names)

            location = jobname + '/graph_frac_sequence.dat'
            destination = base_dir + '/compiled_data/output_' + str(exp_num) + '/'
            if os.path.exists(location):
                shutil.move(location, destination)
                new_name = destination + 'graph_frac_sequence_' + str(ijob) + '.dat'
                old_name = destination + 'graph_frac_sequence.dat'
                os.rename(old_name,new_name)
                with open(inputrecord, 'a') as f:
                    param_names = "Moved" + location + "to " + destination + " and renamed it" + new_name
                    f.write(param_names)

        # Finally, clear out the huge amount of dfn output files.
        if os.path.exists(jobname):
            print("removing all unnecessary dfnWork output files...")
            shutil.rmtree(jobname)

    print(f'made to end of the parallelized function for job {ijob}')
    with open(inputrecord, 'a') as f:
        param_names = f"made it to the end of the parallelized function for job {ijob}"
        f.write(param_names)
    return 0 
   
data = list(range(num_of_experiments*num_of_replications))

max_cpu = 8
num_cpu = min(max_cpu, len(data))
pool = mp.Pool(num_cpu)
outputs = pool.map(run_job, data, chunksize=1)
pool.close()
pool.join()
pool.terminate()



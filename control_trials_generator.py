import random
import os

logo = """
                                  ____
                                 /\\' .\    _____
         C o n t r o l          /: \___\  / .  /\\
          T r i a l s           \\' / . / /____/..\\
       G e n e r a t o r         \/___/  \\'  '\  /
               ( beta )                   \\'__'\/
"""

print(logo)
print()

project_name = str(input(' (✿◠‿◠)  Please enter a name for your experiment:  '))
batch_size = int(input(' (✿◠‿◠)  Please tell me how many MIM numbers would you like to generate per batch:  '))
num_batches  = int(input(' (✿◠‿◠)  Please tell me how many batches you would like to generate:  '))
experimental_MIM_input = str(input(' (✿◠‿◠)  Please list the experimental MIMs to run against controls, separated by spaces (or press RETURN for none):  '))
API_key = str(input(' (✿◠‿◠)  Please enter your API key so I can run the experiment:  '))
print(experimental_MIM_input)
if experimental_MIM_input != '':
    experimental_MIM_list = experimental_MIM_input.split(" ")
else:
    experimental_MIM_list = []

if len(experimental_MIM_list) > 0:
    print(experimental_MIM_list)
#project_name = 'hypophosphatasia'
project_name = project_name.upper()


my_file = open("ALL_OMIM_MIM_NUMBERS_2022-06-18.txt", "r")
MIM_list = my_file.readlines()
my_file.close()

for i in range(len(MIM_list)):
    MIM_list[i] = MIM_list[i].replace("\n","")
    MIM_list[i] = int(MIM_list[i])

if len(experimental_MIM_list) > 0:
    for i in experimental_MIM_list:
        MIM_list.remove(int(i))

print()
print('\(✿◠‿◠)/  I\'m generating batches and running the experiment!')
print()

if os.path.isdir('../'+project_name+'_EXPERIMENT'):
    pass
else:
    os.mkdir('../'+project_name+'_EXPERIMENT')

for i in range(int(num_batches)):
    batch = random.sample(MIM_list, batch_size)
    if len(experimental_MIM_list) > 0:
        for j in experimental_MIM_list:
            batch += [j]
    print(batch)
    f = open('../' + project_name + '_EXPERIMENT/MIM_batch_' + str(i)+'.txt', "w")
    for j in batch:
        f.write(str(j)+'\n')
    f.close()
    for i in range(num_batches):
        f = open('../'+ project_name + '_EXPERIMENT/experiment_'+ str(i) + '.pbs', "w")
        f.write("""#!/bin/bash
#PBS -N """ + project_name + """
#PBS -l select=1:interconnect=fdr:ncpus=20:mem=100gb,walltime=60:00:00

source activate GenoPheno
cd /scratch1/dkeane2/GenoPheno
cp -R ../GenoPheno ../"""+project_name+"""_EXPERIMENT/GenoPheno_""" + str(i) + """
cd ../"""+project_name+"""_EXPERIMENT/GenoPheno_""" + str(i) + """
./GenoPheno.sh ../MIM_batch_""" + str(i)+""".txt """+API_key+""" """ + project_name + """_""" + str(i))
    f.close()

os.chdir('../')


#bashCommand1 = 'git clone http://github.com/devinkeane/GenoPheno'
#os.system(bashCommand1)

for i in range(num_batches):
    bashCommand2 = 'qsub ./'+ project_name +'_EXPERIMENT/experiment_'+ str(i) + '.pbs'
    os.system(bashCommand2)


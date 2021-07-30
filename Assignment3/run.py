import os
import sys

os.system('make clean')
os.system('make')
os.system('bash checkhosts.sh')

file = open("hostfile","r") 
data = file.read() 
hosts = data.split("\n")
if len(hosts)<2:
    os.system('make clean')
    print("Minimum number of hosts required to run the program are not available...exiting now")
    exit()

nodes = [1, 2]
ppns = [1, 2, 4]
for itr in range(5):
    print ("Iteration " + str(itr+1) + " starting...")
    for node in nodes:
        for ppn in ppns:
            rounds = 50 if node is 2 else 200
            string='mpiexec -np '+ str(ppn*node) + ' -ppn '+ str(ppn) + ' -f hostfile ./code tdata.csv '+ str(rounds)
            os.system(string)
            print ("    Command with nodes = " + str(node) + " and ppn = " + str(ppn) + " has executed successfully")
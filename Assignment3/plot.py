import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

sns.set()

intra_node = pd.DataFrame.from_dict({
 	"ppn":[],
    "time":[],
})

inter_node = pd.DataFrame.from_dict({
 	"ppn":[],
    "time":[],
})

f = open(f'dump', "r")
data = f.readlines()
line = 0

for execution in range(5):
    for nodes in [1, 2]:
        for ppn in [1, 2, 4]:
            nums = data[line]
            if nodes is 1:
                intra_node = intra_node.append({
                "ppn": ppn, "time": float(nums)}, ignore_index=True)
            elif nodes is 2:
                inter_node = inter_node.append({
                "ppn": ppn, "time": float(nums)}, ignore_index=True)
            
            line += 1

g = sns.barplot(x = "ppn", y = "time", data = intra_node)
g.set(xlabel='ppn', ylabel='time')
plt.title("Nodes = 1", fontsize=15)
plt.savefig(f'plot_intranode.jpg')

f = sns.barplot(x = "ppn", y = "time", data = inter_node)
f.set(xlabel='ppn', ylabel='time')
plt.title("Nodes = 2", fontsize=15)
plt.savefig(f'plot_internode.jpg')
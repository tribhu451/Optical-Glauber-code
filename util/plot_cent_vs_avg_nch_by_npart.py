import matplotlib.pyplot as plt
import sys


file1 = open(sys.argv[1], 'r')
Lines = file1.readlines()

count = 0
X1, Y1   = [],[]
for line in Lines:
	count += 1
	if count > 1 :
		values = [float(s) for s in line.split()]
		X1.append((values[0]+values[1])/2.)                       
		Y1.append(values[4]) 

	else :
		continue
                    



#ASE1 = [LE1, UE1]                  

file1 = open(sys.argv[2], 'r')
Lines = file1.readlines()

count = 0
X2, Y2   = [],[]
for line in Lines:
	count += 1
	if count > 1 and count < 8 :
		values = [float(s) for s in line.split()]
		X2.append(values[0])                       
		Y2.append(values[8]) 

	else :
		continue
                    

## xlim ylim
fig = plt.figure()
ax = fig.gca()
plt.xlim(-0.5,55.5)
plt.ylim(0.5,9.0)



#plotting ... 


plt.scatter(X1, Y1,marker = '*', s =50, label = r'$  \langle  \frac{dN_{ch}}{d \eta} \rangle | _{\eta = 0} \  /  \ \langle  \frac{dN_{ch}}{d \eta} \rangle _{35-45\%} | _{\eta = 0} $(PHOBOS DATA)', color = 'black',linewidth=0.5)
plt.plot(X2, Y2, ls = ':', label = r'Two component model', color = 'red', lw = 2.5)



#plot label ...
plt.xlabel(r'centrality (%)',fontsize = 15)
plt.ylabel(r'$\langle  \frac{dN_{ch}}{d \eta} \rangle$ ratio',fontsize = 15)
plt.yscale('linear')
plt.legend(fontsize = 12)



#some texts in plot
textstr1 = 'optical glauber \n $x_{hard} = 0.14 $ \n $ds/dy |_{y = 0} \propto $ $(1-X_{hard}) N_{part}/2 + X_{hard}N_{coll}$ '

# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='yellow', alpha=0.5)
ax.text(0.25, 0.75, textstr1, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props, color = 'blue')





ax.tick_params(bottom=False, top=True, left=True, right=True)
#ax.tick_params(labelbottom=False, labeltop=True, labelleft=True, labelright=True)


#plt.show()
plt.savefig('xhard_for_cent_detemination3.pdf', format = 'pdf')
                
    

    


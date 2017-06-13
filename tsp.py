#from random import randint
import numpy as np
import matplotlib.pyplot
import pylab
from math import sqrt
from random import shuffle, randint

def distance(xpt1, xpt2, ypt1, ypt2):
    dist = sqrt((xpt1 - xpt2)**2 + (ypt1 - ypt2)**2)
    return dist

def total_distance(order_list):
    total_dist = 0.
    for i,j in zip(order_list, order_list[1:]):
        total_dist += distance(x[i], x[j],y[i],y[j])
    return total_dist
        
def plot_trip(order_list):
    matplotlib.pyplot.scatter(x,y)
    matplotlib.pyplot.scatter(x[0], y[0], color="red")
    for i, txt in enumerate(city_names):
        matplotlib.pyplot.annotate(txt, (x[i]+2, y[i]+2))
    
    for i,j in zip(order_list, order_list[1:]):
        matplotlib.pyplot.plot([x[i], x[j]],[y[i],y[j]])
    #matplotlib.pyplot.plot(x,y)
    matplotlib.pyplot.show()        

def crossover(order1, order2):
    order_length = len(order1)
    start_int = randint(0,order_length/2-1)
    new_order = order1[start_int:start_int+order_length/2]
    append_list = list(set(order2) - set(new_order)) #sets are unordered by nature
    append_list_ordered = [o for o in order2 if o in append_list] #so rearrange to original order
    new_order[order_length/2:-1] = append_list_ordered #add genetic from order2 to offspring
    return new_order

n_cities = 10
#x = [np.random.randint(1,100) for number in range(n_cities)]
#y = [np.random.randint(1,100) for number in range(n_cities)]

x = [33, 44, 56, 1, 34, 67, 98, 64, 22, 71]
y = [10, 71, 22, 44, 23, 99, 56, 67, 42, 34]
city_names = [0,1,2,3,4,5,6,7,8,9]
order = [0,1,2,3,4,5,6,7,8,9]
order_next = [4,5,6,3,2,7,8,1,9,0]
crossover(order,order_next)
#best_order = [0,4,3,8,1,7,5,6,9,2]

current_best = float("inf")



#plot_trip(best_order)
#print total_distance(best_order)

#best_distance = total_distance(order)
#for i in range(2):
#    shuffle(order)
#    if total_distance(order) < best_distance:
#        best_order = order
#    print i

#print total_distance(best_order)
#plot_trip(best_order)

from functools import lru_cache
#from backports.functools_lru_cache import lru_cache
import matplotlib.pyplot as plt 

plt.style.use('dark_background')

@lru_cache(maxsize = 50000)
def fibonacci(n):
    if n == 1:
        return 1
    elif n<0 :
        raise ValueError("Value must be Sup or equal 1")
    elif n==2:
        return 1
    elif n>2:
        return fibonacci(n-1)+fibonacci(n-2)

if __name__=='__main__':
    fig =  plt.figure()
    ax=fig.add_subplot(1,1,1)
    x_Data = list(range(1,101))
    y_Data = []
    for i in x_Data:
        y_Data.append(fibonacci(i))
    
    ax.plot(x_Data,y_Data,'w',linestyle='dotted',marker='o',markerfacecolor='red',markersize=7)
    plt.show()
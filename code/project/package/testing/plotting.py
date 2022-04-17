from ..helpers import *

def plot_results(x, y, plot_name = None, x_axis_name = None, y_axis_name = "Average State Fidelity"):
    
    plt.scatter(x, y)
    
    if plot_name:
        plt.title(plot_name)
        
    if x_axis_name:
        plt.xlabel(x_axis_name)
    if y_axis_name:
        plt.ylabel(y_axis_name)
        
    plt.show()



	
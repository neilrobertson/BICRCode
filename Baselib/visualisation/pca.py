'''
Created on 23 Aug 2010

@author: mcbryan
'''
from numpy import *
import pylab
import csv
#from rpy import *
from rpy2.rpy_classic import *
from matplotlib.pyplot import *

def pcaPlot(col_names, array_colour, outputFolder, componentsfile = None):
    
    if componentsfile == None:
        componentsfile = outputFolder+"/principalcomponents.csv"
    
    with open(componentsfile) as infile:
        components = csv.reader(infile, delimiter="\t", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
        pca_names = None
        pcas = {}
        for row in components:
            if pca_names == None:
                pca_names = row
            else:
                pcas[row[0]] = row[1:]
        
        
        figure(1)
        for i in range(len(pca_names)):
            
            xpadding = 0.1
            xpos = 0.0+xpadding
            width = 0.7
            
            rows = len(pca_names) + 1
            
            ypadding = 0.05
            height = (1.0-ypadding)/rows - ypadding/2
            ypos = 1.0 - ((float(i)-ypadding)/rows + ypadding/2) - height
            
            ax = axes([xpos, ypos, width, height])
            # plot
            component = []
            
            # this component id
            for col_name in col_names:
                component.append(pcas[col_name][i])
            stem(arange(0, len(col_names)), component, '-.',  
                 #markersize=20
                 )
            
            axhline(y=0)
            xlim(-0.5, len(col_names)-0.5)
            yticks(visible=True, size=4)
            ylabel(pca_names[i])
            
            if i == len(col_names)-1:
                xticks(arange(0, len(col_names)), col_names, visible=True,  size=5, rotation=15)
            else:
                xticks(visible=False)
            
            ylims = ylim()
            
            ax = axes([xpos+width, ypos, 0.1, height])
            for col_name in col_names:
                plot(pcas[col_name][i], '.',  markersize=20, color=array_colour(col_name), label=col_name)
                xlim(0, 0)
                xticks(visible=False)
                yticks(visible=False)
                ylim(ylims)
        
        #handles, labels = ax.get_legend_handles_labels()
        #figlegend(handles,  labels, "lower center", mode="expand", ncol=4)
        
        savefig( open(outputFolder+"/pca-1d.png", "w"), format='png' )
        pylab.clf()
        
        
        fig = figure(figsize=(16,12))
        
        fig.subplots_adjust(bottom=0.2, top=0.95, left=0.05, right=0.95)
        
        markers =  [        'o', #  circle
                            's', #  square
                            '^', #  triangle up
                            '>', #  triangle right
                            'v', #  triangle down
                            '<', #  triangle left
                            'd', #  diamond
                            'p', #  pentagram
                            'h', #  hexagon
                            
                            ] * 10
        
        ax = axes()
        
        for col_name in col_names:
            scatter(pcas[col_name][0], pcas[col_name][1], s=50, c=array_colour(col_name), marker=markers[col_names.index(col_name)], label=col_name)
        
        xlabel(pca_names[0])
        ylabel(pca_names[1])
        
        handles, labels = ax.get_legend_handles_labels()
        
        figlegend(handles,  labels, "lower center", mode="expand", ncol=4, prop={"size":10})
        
        fig.savefig( open(outputFolder+"/pca-2d.png", "w"), format='png')
        pylab.clf()
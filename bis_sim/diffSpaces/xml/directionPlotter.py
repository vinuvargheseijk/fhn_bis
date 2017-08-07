import sys
from matplotlib import cm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import xml.etree.ElementTree as ET
from scipy import stats
import glob
import os
import re

if len(sys.argv) < 2:
    print "Usage: ", sys.argv[0], " outputfilename" " searchTerm" " [prefix to plotTitle]"
    quit() 
outputFileName = str(sys.argv[1])
searchTerm = str(sys.argv[2])

xmlFileList = []
for file in glob.glob("*.xml"):
    xmlFileList.append(file)	# to get all the xml files in the directory
    print file
n = len(xmlFileList) 	# n will determine how many subplots the plot is to be divided into



def plotDirectionMetric( fig, tree, index, plotTitle ):
    # Note that we are transposing X and Y axes for the matrix heatmap plot.
    #showYlabel = (index == 0) or (index == 3) or (index == 6)
    xaxis = tree.find( 'xaxis' )
    assert( xaxis != None )
    xtitle = xaxis.get( 'title' )
    
    yaxis = tree.find( 'yaxis' )
    assert( yaxis != None )
    ytitle = yaxis.get( 'title' )

    xValues = [float(j) for j in yaxis.text.split()]
    yValues = [float(j) for j in xaxis.text.split()]
    
    seqOutput = tree.find( 'seqOutput' )
    stimSequence = tree.findall( 'seqOutput/stimSequence' )
    score = { int(i.get( 'index' )):float(i.get('score')) for i in stimSequence}
    
    orderedScore = [ score[i] for i in range( len( score ) ) ]
    labelIndex = 10

    
    #################### color heatmap here ###################
    simOutput = tree.find( "simOutput" )
    aocSumAll = []
    for dtData in simOutput.iter( "dtData" ):
        for distanceData in dtData.iter( "distanceData" ):
            distance = int( distanceData.get( 'distance' ) )
            for labelData in distanceData.iter( "labelData" ):
                label = labelData.get( 'title' )
                if label == 'aocSum':
                    values = [float(q) for q in labelData.text.split()]	# [values] stores the aocSum values
                    # aocSumAll.append( (values[0] - np.mean( values) ) / max( values ) )
                    aocSumAll.append( (values[0] - values[-1]) / max( values[0],values[-1] ) )
    ax = plt.subplot(3,3,index+1)
    ax.tick_params( direction = 'out' )
    d = np.array( aocSumAll)
    md = min(0.0, min( d ))
    # This has two problems: The axes are transposed, and the ordering of
    # the entries for time (in the columns; xValues) are flipped.
    d.shape = (len( xValues ), len( yValues ) )
    td0 = [ i[::-1] for i in d ] # Unflip the time ordering
    td = np.transpose( td0 )     # Put time on the x axis, dist on y axis.
    cax = plt.imshow( td, cmap = cm.cool, vmin = md, interpolation='none' )
    cbar = fig.colorbar( cax, fraction=0.046, pad=0.04 )
    if len(sys.argv) == 4:
    	plt.title(str(sys.argv[3] + plotTitle))
    else:
    	plt.title(plotTitle)
    plt.xlabel( "Time (s)", fontsize = 14 )
    plt.ylabel( r"Distance ($\mu$m)", fontsize = 14)
    xticks = [str(int(x)) for x in xValues]
    plt.xticks ( range( len(xValues ) ), xticks )
    plt.yticks( range( len(yValues)-1, -1, -1 ), [str(y) for y in yValues] )

def main():
    fig = plt.figure( figsize = (20.5,10.5), facecolor='blue')
    for index,xmlFile in enumerate(xmlFileList):
    	bisTree = ET.parse( xmlFile )
    	plotTitle = str(xmlFile).replace(searchTerm,'').replace('.xml','').replace('_','-')
    	plotDirectionMetric( fig, bisTree, index, plotTitle )
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    # fig.subplots(squeeze=False)
    plt.suptitle( outputFileName + ".png", fontsize=20 )

    plt.savefig( outputFileName + ".png" )

if __name__ == '__main__':
    main()

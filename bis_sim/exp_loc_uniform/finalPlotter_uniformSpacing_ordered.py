import sys
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import xml.etree.ElementTree as ET
from scipy import stats

if len(sys.argv) < 2:
    print "Usage: ", sys.argv[0], " xmlfilename"
    quit() 
outXmlFile = str(sys.argv[1])
dt = str(sys.argv[2])


def plotScatter( fig, tree, dist, dt, index ):
    # Note that we are transposing X and Y axes for the matrix heatmap plot.
    showXlabel = (index == 3)
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
    
    scatter = tree.findall("./simOutput/dtData[@dt='" + dt + "']/distanceData[@distance='" + dist + "']")
    labelIndex = 10

    #panelName = chr( ord( 'A' ) + index )
    panelName = 'd = ' + str(dist) + ',\ndx = ' + str(int(dist)/5) + ',\ndt = ' + str(dt)
    print "Lengths: scatter = ", len(scatter), 
    if len( scatter ) == 1:
        labelData = scatter[0].findall( 'labelData' )
        print ", labels = ", len( labelData )
        if len( labelData ) > 0:
            #plt.suptitle( 'coords = (' + dist + ', ' + dt + ')', fontweight = 'bold' )
            ax = plt.subplot( 4, 2, 2 * index + 1 )
            data = np.array( [float(j) for j in labelData[labelIndex].text.split()] ) * 0.02 # 0.02 is the chemPlotDt from the simulations.
            assert( len(data) == len(orderedScore))
            slope, intercept, r_value, p_value, std_err = stats.linregress( orderedScore, data )
            cax = plt.scatter( orderedScore, data )
            ax.spines['top'].set_visible( False )
            ax.spines['right'].set_visible( False )
            for tick in ax.xaxis.get_major_ticks():
                tick.tick2On = False
            for tick in ax.yaxis.get_major_ticks():
                tick.tick2On = False
            ax.tick_params( direction = 'out' )
            ax.xaxis.set_ticks( (-1, 0, 1) )
            ax.text( -0.32, 1, panelName, fontsize = 18, weight = 'bold', transform=ax.transAxes )
            #slopeLabel = "m={0:.2f}, R^2={1:.2f}".format( slope, r_value**2)
            #plt.text( min( orderedScore ), max( data ), slopeLabel )
            #plt.title( labelData[labelIndex].get( 'title' ) )
            if showXlabel:
                plt.xlabel( 'Sequence Score', fontsize = 14 )
            plt.ylabel( 'Total A', fontsize = 14 )

    #################### color heatmap here ###################
    simOutput = tree.find( "simOutput" )
    aocSumAll = []
    for dtData in simOutput.iter( "dtData" ):
        for distanceData in dtData.iter( "distanceData" ):
            distance = int( distanceData.get( 'distance' ) )
            for labelData in distanceData.iter( "labelData" ):
                label = labelData.get( 'title' )
                if label == 'aocSum':
                    values = [float(q) for q in labelData.text.split()]
                    aocSumAll.append( (values[0] - np.mean( values) ) / max( values ) )
    ax = plt.subplot( 4,2,2*0 + 2 )
    ax.tick_params( direction = 'out' )
    d = np.array( aocSumAll)
    md = min( 0.0, min(d) )
    # This has two problems: The axes are transposed, and the ordering of
    # the entries for time (in the columns; xValues) are flipped.
    d.shape = (len( xValues ), len( yValues ) )
    td0 = [ i[::-1] for i in d ] # Unflip the time ordering
    td = np.transpose( td0 )     # Put time on the x axis, dist on y axis.
    cax = plt.imshow( td, cmap = cm.cool, vmin = md, interpolation='none' )
    cbar = fig.colorbar( cax )
    #plt.title ( "The hippos are coming out to play" )
    if showXlabel:
        plt.xlabel( "Time (s)", fontsize = 14 )
    plt.ylabel( r"Distance ($\mu$m)", fontsize = 14)
    xticks = [str(int(x)) for x in xValues]
    # xticks[0] = '0.5'
    # xticks[2] = '1.5'
    #plt.xticks ( range( len(xValues ) ), [str(int(x)) for x in xValues] )
    plt.xticks ( range( len(xValues ) ), xticks )
    plt.yticks( range( len(yValues)-1, -1, -1 ), [str(y) for y in yValues] )
    #fig.tight_layout()




def main():
    #fig = plt.figure(figsize = (7,10), facecolor='white')
    fig = plt.figure( figsize = (18.5,10.5), facecolor='blue')
    bisTree = ET.parse( outXmlFile )
    #bisTree = ET.parse( 'bis4.0.xml' )
    #fhnTree = ET.parse( 'fhn7.0.xml' )
    #negFBTree = ET.parse( 'negFB6.1.xml' )
    #negFFTree = ET.parse( 'negFF9.1.xml' )

    plotScatter( fig, bisTree, '10', dt, 0 )
    plotScatter( fig, bisTree, '15', dt, 1 )
    plotScatter( fig, bisTree, '20', dt, 2 )
    plotScatter( fig, bisTree, '25', dt, 3 )
    #plt.show()
    #ampl = outXmlFile[-8:].replace('.xml','').replace('_','')
    plt.savefig("fhn_uniformSpacing_" + "dt-" + dt + '_tuningOrdered' + '.png')

if __name__ == '__main__':
    main()

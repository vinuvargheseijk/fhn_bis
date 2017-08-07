import moose
import pylab
import rdesigneur as rd
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
import rdesigneur as rd
import xml.etree.ElementTree as ET
import itertools
from scipy import stats

params = {
    'diffusionLength':0.5e-6,  # Diffusion characteristic length, used as voxel length too.
    'dendDiameter': 1.0e-6, # Diameter of section of dendrite in model
    'dendLength': 0.5e-6,    # Length of section of dendrite in model
    'diffConstCa':100e-12,  # Diffusion constant of Ca
    'stimAmplitude': 0.005, # Ca Stimulus amplitude, mM
    'baseCa':2.5e-4,        # Base Ca level, mM.
    'BAPCa':0.002,          # Dend Ca spike amplitude
    'BAPwidth':0.1,         # Dend Ca spike width.
    'blankVoxelsAtEnd':10,  # of voxels to leave blank at end of cylinder
    'preStimTime':10.0,     # Time to run before turning on stimulus.
    'postStimTime':40.0,    # Time to run after stimulus.
    'stimWidth': 2.9,       # Duration of Ca influx for each stimulus.
    'spineSpacing':1.1e-6,  # Spacing between spines.
    'diffConstMAPK': 5e-12, # Diffusion constant for MAPK
    'diffConstPP': 2e-12,   # Diff constant for MAPK-activated phosphatase
    'CaActivateRafKf': 6e6, # 1/sec/mM^2: rate for activation of Raf by Ca
    'cellModel':'PassiveSoma',  # Cell morphology script
    'chemModel':'NN_mapk14.g',  # Chem model definition
    #'seqDt': 3.0,           # Time interval between successive inputs in seq
    'seqDt' : 3.0,
    'seqDx': 3.0e-6,        # Distance between successive inputs in seq.
    #'seed': 3,            # Seed for random number generator
    'seed': 12345,           # Seed for random number generator
    'sequence': '01234',    # Sequence of spines, spaced by seqDx microns, 
                            # activated every seqDt seconds
    'fnumber': 0,           # identifier for run
}

RM=1.0
RA=10.0
CM=0.001
lenbfNotch=3
lenafNotch=3
displayMoogli = False

def setDiffConst( element, paramName ):
    e = moose.element( '/library/chem/dend/DEND/' + element )
    e.diffConst = params[ paramName ] 
    
def runStimulus( sequence ):
    baseCa = params['baseCa'] / 50.0
    stimQ = buildStimulusQ( baseCa, sequence )
    Ca_input = moose.vec( '/model/chem/psd/Ca_input' )
    Ca_input.concInit = baseCa
    moose.vec( '/model/chem/dend/DEND/Ca_input' ).concInit = baseCa
    moose.reinit()
    
    clock = moose.element( '/clock' )
    for t in sorted( stimQ ):
        [index, conc] = stimQ[t]
        currt = clock.currentTime
        if ( t > currt ):
            moose.start( t - currt )
            print "At t = ", t, "; assigning CaInput[", index, "] = ", conc
            Ca_input[ index ].concInit = conc
    moose.start(params['postStimTime'] )
    print "Finished stimulus run at t = ", clock.currentTime

def buildStimulusQ( baseCa, sequence ):
    #blanks = params['blankVoxelsAtEnd']
    #step = int( round( params['seqDx'] / params['spineSpacing'] ) )
    sequence = [ int(i) for i in sequence ]
    stimulusQ = {}
    onCa = params['stimAmplitude']
    stimStart = params['preStimTime']
    for i in sequence:
        #stimulusQ[ stimStart ] = [ blanks + i * step, onCa ]
        stimulusQ[ stimStart ] = [i, onCa ]
        stimEnd = stimStart + params['stimWidth']
        #stimulusQ[ stimEnd ] = [ blanks + i * step, baseCa ]
        stimulusQ[ stimEnd ] = [i, baseCa ]
        stimStart += params['seqDt']
    return stimulusQ

def runAndDisplay( seq, label, pos, chemPlotDt, maxy ):
    runStimulus( seq )
    mapk = moose.vec( '/model/graphs/plot0' )
    parentVoxel = moose.element( '/model/chem/spine' ).neuronVoxel
    blanks = params['blankVoxelsAtEnd']
    #step = int( round( params['seqDx'] / params['spineSpacing'] ) )
    #plotIndices = [int(parentVoxel[blanks + i * step]) for i in range(len(seq)) ]
    plotIndices = [int(parentVoxel[i]) for i in range(len(seq)) ]
    print plotIndices
    t = np.arange( 0, len( mapk[0].vector ) * chemPlotDt, chemPlotDt )
    ax = plotBoilerplate( label, pos, 'Time (s)', 'MAPK-P ($\mu$M)')
    maxx = max( t )
    xt = np.arange( 0, maxx, 20 )
    xticks = [ str(int(i)) for i in xt ]
    print "XT = ", xt, xticks

    ax.xaxis.set_ticks( xt )
    ax.set_xticklabels( xticks )
    ax.set_ylim( [0.0, maxy] )
    print 'size of MAPK is ', len(mapk[1].vector), plotIndices, 'size of t is', len(t)
    print 'size of MAPK i is ', len(mapk)
    #for i in range(0,len(plotIndices),1):
    for i in plotIndices:
        print 'i is ', i
        print mapk[i].vector
        plt.plot( t, mapk[i].vector * 10000 )

def panelEFspatialSeq( fig ):
    print "Starting Panel EF"
    moose.seed( int(params['seed']) )
    rdes = rd.rdesigneur(
        useGssa = False,
        turnOffElec = True,
        chemPlotDt = 0.02,
        diffusionLength = params['diffusionLength'],
       # spineProto = [['makePassiveSpine()', 'spine']],
       # spineDistrib = [['spine', '#', str(params['spineSpacing']),'1e-7','1.4','0']],
        cellProto = [['./non_branch.p', 'elec']],
        chemProto = [[params['chemModel'], 'chem']],
        chemDistrib = [['chem', '#', 'install', '1' ]],
        plotList = [
            ['#', '1', 'dend/DEND/P_MAPK', 'conc', '[dend P_MAPK]'],
            ['#head#', '1', 'psd/Ca', 'conc', '[PSD Ca]'],
            ['#head#', '1', 'spine/Ca', 'conc', '[spine Ca]'],
            ['soma', '1', 'dend/DEND/Ca', 'conc', '[dend Ca]'],
            ],
    )   
    # Assign parameters to the prototype model.
    setDiffConst( 'Ca', 'diffConstCa' )
    setDiffConst( '../../spine/Ca', 'diffConstCa' )
    setDiffConst( '../../psd/Ca', 'diffConstCa' )
    setDiffConst( 'P_MAPK', 'diffConstMAPK' )
    setDiffConst( 'MAPK', 'diffConstMAPK' )
    setDiffConst( 'reg_phosphatase', 'diffConstPP' )
    setDiffConst( 'inact_phosphatase', 'diffConstPP' )
    moose.element( '/library/chem/dend/DEND/Ca_activate_Raf' ).Kf = params['CaActivateRafKf']
    print "Set up rdesigneur"

    rdes.buildModel()
    print "MODEL BUILT"

    ################################################################
    # Run and display the stimulus
    runAndDisplay( '01234', 'E', (1,0), rdes.chemPlotDt, 7 )
    runAndDisplay( '40312', 'F', (1,1), rdes.chemPlotDt, 7 )

    #moose.delete( '/model' )

def plotBoilerplate( panelTitle, plotPos, xlabel = '', ylabel = '', xticks = [] ):
    ax = plt.subplot2grid( (5,2), plotPos )
    #ax.xaxis.set_ticks( i[1] )
    #ax.locator_params( 
    ax.spines['top'].set_visible( False )
    ax.spines['right'].set_visible( False )
    ax.tick_params( direction = 'out' )
    #ax.set_xticklabels([])
    ax.set_xticklabels( xticks )
    ax.set_xlabel( xlabel )
    for tick in ax.xaxis.get_major_ticks():
        tick.tick2On = False
    for tick in ax.yaxis.get_major_ticks():
        tick.tick2On = False

    plt.ylabel( ylabel, fontsize = 14 )
    # alternate way of doing this separately.
    #plt.yaxis.label.size_size(16)
    #plt.title( 'B' )
    ax.text( -0.3, 1.1, panelTitle, fontsize = 18, weight = 'bold',
    transform=ax.transAxes )
    return ax


    
def main():
    global params
    fig = plt.figure(figsize = (10,6), facecolor='white')

    library = moose.Neutral( '/library' )

    for ii in range( len( sys.argv ) ):
        if sys.argv[ii][:2] == '--':
            argName = sys.argv[ii][2:]
            if argName in params:
                params[argName] = float( sys.argv[ii+1] )
                if argName == 'sequence':
                    params[argName] = sys.argv[ii+1] # Leave it as a str.

    moose.seed( int(params['seed']) )
    panelEFspatialSeq( fig );
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()


    


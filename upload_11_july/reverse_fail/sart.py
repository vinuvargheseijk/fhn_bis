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
    'tau':2.5,  # Time-const.
    'a':0.7,  # Coeff1
    'b':0.8,  # Coeff2
    'offsetV': 2.0, # Offset for V for FHN eqns.
    'offsetW': 0.8, # Offset for W for FHN eqns.
    'diffusionLength':1.0e-6,  # Diffusion characteristic length, used as voxel length too.
    'dendDiameter': 10e-6,  # Diameter of section of dendrite in model
    'dendLength': 100e-6,   # Length of section of dendrite in model
    'diffConstA':9.0,      # Diffusion constant of A
    'diffConstB':5.0,       # Diffusion constant of B
    'stimWidth' :1.0,        # Stimulus width in seconds
    'stimAmplitude':0.4,    # Stimulus amplitude, arb units. From FHN review
    'blankVoxelsAtEnd':20,  # of voxels to leave blank at end of cylinder
    'preStimTime':10.0,     # Time to run before turning on stimulus.
    'postStimTime':40.0,    # Time to run after stimulus. ~3x decay time
    'settleTime':20.0,    # Settling time of response, after stimulus. 
                          # To include in analysis of total response over 
                          # entire dendrite.
    'fnumber': 1,         # Number to append to fname
}

RM=1.0
RA=10.0
CM=0.001
lenbfNotch=3
lenafNotch=3
displayMoogli = False

def sp( arg, term ):
    return str( params[arg] ) + term

def setDiffConst( element, paramName ):
    e = moose.element( '/library/chem/dend/DEND/' + element )
    e.diffConst = params[ paramName ] 
    
def makeChemProto( name, stimAmpl = 1, diffLength = 1e-6, preStim = 10.0 ):
    # Parameters

    sw = params['stimWidth']
    print diffLength
    dca = params['diffConstA'] * diffLength * diffLength
    dcb = params['diffConstB'] * diffLength * diffLength

    # Objects
    chem = moose.Neutral( '/library/' + name )
    compt = moose.CubeMesh( chem.path + '/dend' )
    A = moose.Pool( compt.path + '/A' )
    B = moose.Pool( compt.path + '/B' )
    Z = moose.BufPool( compt.path + '/Z' )
    Ca = moose.BufPool( compt.path + '/Ca' )
    phase = moose.BufPool( compt.path + '/phase' )
    vel = moose.BufPool( compt.path + '/vel' )
    ampl = moose.BufPool( compt.path + '/ampl' )
    Adot = moose.Function( A.path + '/Adot' )
    Bdot = moose.Function( B.path + '/Bdot' )
    CaStim = moose.Function( Ca.path + '/CaStim' )
    A.diffConst = dca
    B.diffConst = dcb

    # Equations

    xv0 = '(x0-' + str( params['offsetV'] ) + ')'
    xv1 = '(x1-' + str( params['offsetV'] ) + ')'
    xw1 = '(x1-' + str( params['offsetW'] ) + ')'
    xw2 = '(x2-' + str( params['offsetW'] ) + ')'
    Adot.expr = '5*x3*( 0.0 +' + xv1 + ' - (' + xv1 + '^3)/3 -' + xw2 + ' + x0 )' # x0=Ca, x1=A, x2=B
    #Adot.expr = 'x3*( 0.5 + x1 - x1*x1*x1/3 - x2 + x0 )' # x0=Ca, x1=A, x2=B
    #Bdot.expr = 'x2*' + str(1.0/params['tau']) + '*(x0 + ' + sp('a', ' - ') + sp( 'b', ' * x1 ' ) + ')'
    Bdot.expr = 'x2*' + str(1.0/params['tau']) + '*(' + xv0 + '+' + sp('a', ' - ') + sp( 'b', ' * ' + xw1 ) + ')'
    CaStim.expr = 'x2 * exp( -((x0 - t)^2)/(2* ' + str(sw*sw) + ') )'

    print Adot.expr
    print Bdot.expr
    print CaStim.expr

    # Connections
    Adot.x.num = 4
    moose.connect( Ca, 'nOut', Adot.x[0], 'input' )
    moose.connect( A, 'nOut', Adot.x[1], 'input' )
    moose.connect( B, 'nOut', Adot.x[2], 'input' )
    moose.connect( Z, 'nOut', Adot.x[3], 'input' )
    moose.connect( Adot, 'valueOut', A, 'increment' )

    Bdot.x.num = 3
    moose.connect( A, 'nOut', Bdot.x[0], 'input' )
    moose.connect( B, 'nOut', Bdot.x[1], 'input' )
    moose.connect( Z, 'nOut', Bdot.x[2], 'input' )
    moose.connect( Bdot, 'valueOut', B, 'increment' )

    CaStim.x.num = 3
    moose.connect( phase, 'nOut', CaStim.x[0], 'input' )
    moose.connect( vel, 'nOut', CaStim.x[1], 'input' )
    moose.connect( ampl, 'nOut', CaStim.x[2], 'input' )
    moose.connect( CaStim, 'valueOut', Ca, 'setN' )

    return compt

def runStimulus( sequence, diffusionLength, v ):
    preStim=params['preStimTime']
    baseCa = params['baseCa'] / 50.0
    stimQ = buildStimulusQ( baseCa, sequence )
    vel = moose.vec('/model/chem/dend/vel')
    vel.nInit = v * diffusionLength
    #Ca_input = moose.vec( '/model/chem/psd/Ca_input' )
    #Ca_input.concInit = baseCa
    #moose.vec( '/model/chem/dend/DEND/Ca_input' ).concInit = baseCa
    Z=moose.vec('/model/chem/dend/Z')
    phase = moose.vec('/model/chem/dend/phase')
    ampl = moose.vec('/model/chem/dend/phase')
    ampl.nInit = params['stimAmplitude']
    Z.nInit = 0
    phase.nInit = 1000.0
    moose.reinit()
    
    clock = moose.element( '/clock' )
#    for t in sorted( stimQ ):
#        [index, conc] = stimQ[t]
#        currt = clock.currentTime
#        if ( t > currt ):
#            moose.start( t - currt )
#            print "At t = ", t, "; assigning CaInput[", index, "] = ", conc
#            Z[ index ].nInit = 1
    pos = (11,13,15,17,19)
    temp=pos[1]
    sequence = [int(j) for j in sequence]
    print sequence
    for j in sequence:
       k = pos[j]
       Z[k].nInit=1
       seqDt = params['seqDt']
       #print sequence[j]
       phase[k].nInit=preStim+j*seqDt
       print 'k is', k, 'phase is', phase[k].n
    #print 'phase is', phase[11].n,phase[12].n,phase[13].n,phase[14].n
    #Z[temp].nInit=0
    #phase[11].nInit=10
    #phase[13].nInit=20
    ##phase[15].nInit=30
    moose.reinit()
    runTime=100
    moose.start(runTime)
    
    #moose.start(params['postStimTime'] )
    print "Finished stimulus run at t = ", clock.currentTime

def buildStimulusQ( baseCa, sequence ):
    #blanks = params['blankVoxelsAtEnd']
    #step = int( round( params['seqDx'] / params['spineSpacing'] ) )
    preStim = 10.0
    sequence = [ int(i) for i in sequence ]
    stimulusQ = {}
    onCa = params['stimAmplitude']
    stimStart = params['preStimTime']
    pos = (11,13,15,17,19)
    for i in sequence:
        #print 'i in sequence is', i
        #k = pos[i]
        #phase[k].nInit=preStim+i*seqDt
        #stimulusQ[ stimStart ] = [ blanks + i * step, onCa ]
        stimulusQ[ stimStart ] = [pos[i], onCa ]
        stimEnd = stimStart + params['stimWidth']
        #stimulusQ[ stimEnd ] = [ blanks + i * step, baseCa ]
        stimulusQ[ stimEnd ] = [pos[i], baseCa ]
        stimStart += params['seqDt']
    return stimulusQ

def runAndDisplay( seq, label, pos, chemPlotDt, maxy ):
    diffusionLength = params['diffusionLength']
    v = 3.0
    runStimulus( seq, diffusionLength, v )
    A = moose.vec( '/model/graphs/plot0' )
    
    Ca = moose.vec( '/model/graphs/plot3' )
    print 'size of Ca is', len(Ca)
    phase = moose.vec('/model/graphs/plot2')
   # parentVoxel = moose.element( '/model/chem/spine' ).neuronVoxel
    parentVoxel = [11,13,15,17,19]
    blanks = params['blankVoxelsAtEnd']
    #print size(parentVoxel)
    #step = int( round( params['seqDx'] / params['spineSpacing'] ) )
    #plotIndices = [int(parentVoxel[blanks + i * step]) for i in range(len(seq)) ]
    plotIndices = [int(parentVoxel[i]) for i in range(len(seq)) ]
    
    print plotIndices
    t = np.arange( 0, len( Ca[0].vector ) * chemPlotDt, chemPlotDt )
    
    #t = np.arange( 0, len( phase[0].vector ) * chemPlotDt, chemPlotDt )
    ax = plotBoilerplate( label, pos, 'Time (s)', 'Ca (number of molecules)')
    maxx = max( t )
    xt = np.arange( 0, maxx, 20 )
    xticks = [ str(int(i)) for i in xt ]
    print "XT = ", xt, xticks

    ax.xaxis.set_ticks( xt )
    ax.set_xticklabels( xticks )
    ax.set_ylim( [0.0, maxy] )
    print 'size of A is ', len(A[1].vector), plotIndices, 'size of t is', len(t)
    print 'size of A i is ', len(A)
    #for i in range(0,len(plotIndices),1):
    for i in plotIndices:
        print 'i is ', i
        print 'Ca is', Ca[i].vector
        plt.plot( t, Ca[i].vector)

def panelEFspatialSeq( fig ):
    print "Starting Panel EF"
    moose.seed( int(params['seed']) )
    diffusionLength = params['diffusionLength']
    makeChemProto( 'spread', stimAmpl = 1,diffLength = diffusionLength )
    rdes = rd.rdesigneur(
        useGssa = False,
        turnOffElec = True,
        chemPlotDt = 0.02,
       # spineProto = [['makePassiveSpine()', 'spine']],
       # spineDistrib = [['spine', '#', str(params['spineSpacing']),'1e-7','1.4','0']],
        cellProto = [['./taper.p', 'elec']],
        chemProto = [['spread', 'chem']],
        chemDistrib = [['chem', '#', 'install', '1' ]],
        plotList = [['#', '1', 'dend/A', 'n', '# of A'],
            ['#', '1', 'dend/B', 'n', '# of B'],
            ['#', '1', 'dend/phase', 'n', '# of Phase'],
            ['#', '1', 'dend/Ca', 'n', '# of Ca']
            ],
        moogList = [['#', 'x>5e-6 && x<100e-6', '.', 'Vm', 'Vm']]
        #plotList = [
         #   ['#', '1', 'dend/P_MAPK', 'conc', '[dend P_MAPK]'],
          #  ['#head#', '1', 'psd/Ca', 'conc', '[PSD Ca]'],
           # ['#head#', '1', 'spine/Ca', 'conc', '[spine Ca]'],
            #['soma', '1', 'dend/DEND/Ca', 'conc', '[dend Ca]'],
            #],
    )   
    # Assign parameters to the prototype model.
    #setDiffConst( 'Ca', 'diffConstCa' )
    #setDiffConst( '../../spine/Ca', 'diffConstCa' )
    #setDiffConst( '../../psd/Ca', 'diffConstCa' )
    #setDiffConst( 'P_MAPK', 'diffConstMAPK' )
    #setDiffConst( 'MAPK', 'diffConstMAPK' )
    #setDiffConst( 'reg_phosphatase', 'diffConstPP' )
    #setDiffConst( 'inact_phosphatase', 'diffConstPP' )
    #moose.element( '/library/chem/dend/DEND/Ca_activate_Raf' ).Kf = params['CaActivateRafKf']
    print "Set up rdesigneur"

    rdes.buildModel()
    print "MODEL BUILT"

    ################################################################
    # Run and display the stimulus
    runAndDisplay( '01234', 'E', 1, rdes.chemPlotDt, 1.5 )
    runAndDisplay( '43210', 'F', 2, rdes.chemPlotDt, 1.5 )

    #moose.delete( '/model' )

def plotBoilerplate( panelTitle, plotPos, xlabel = '', ylabel = '', xticks = [] ):
    #ax = plt.subplot2grid( (1,1), plotPos )
    ax = plt.subplot(1, 2, plotPos ) 
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


    


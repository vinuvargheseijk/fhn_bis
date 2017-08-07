import moose
import pylab
import numpy as np
# import numpy
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
import rdesigneur as rd
import xml.etree.ElementTree as ET
import itertools
from scipy import stats

#fname = 'fhn_long_taperModel'
fname = str(sys.argv[2])
displayMoogli = False
moogliDistance = 10
moogliDt = 1.0
mootliSequence = '01234'
displayScatter = False
scatterParamToUseInStats = 0

# Stim amplitude is unitless
# Stim Width is unitless, defined as multiple of diffusion length.
# Stim Vel is unitless, defined in terms of diffusion length by  time units of diff constt.
# diffConst is defined in terms of square of diffusion length by unit time.
# diffLength here is in SI units: m^2/s
#
# Equations here are:
# vdot = v - v^3/3 - w + Iext
# tau.wdot = v + a - bw
# I use A = v and B = w
# I have diffusion only for v (A), as per the original diffusion version.
# Problem with original FHN is that we get negative values.
#

params = {
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
numSpine = 5

def sp( arg, term ):
    return str( params[arg] ) + term

def makeChemProto( name, stimAmpl = 1, diffLength = 1e-6, preStim = 10.0 ):
    # Parameters

    sw = params['stimWidth']
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

def makePassiveSoma( name, length, diameter):
    print("length: ", length)
    print("diameter: ", diameter)
    numComp = 10
    x = length/float(numComp)
    # x = length
    elecid = moose.Neuron( '/library/' + name )
    soma = moose.Compartment( elecid.path + '/soma' )
    soma.diameter = diameter
    soma.length = length
    soma.x0 = 0
    soma.x = x

    ls = x
    n = 2.0
    dendList = []
    for i in range(numComp -1):
        dendList.append(moose.Compartment( elecid.path + '/dend' + str(i) ))
        dend = dendList[-1]
        dend.diameter = diameter/n
        # dl = ls + x
        dend.x0 = ls
        dend.x = ls + x
        dend.length = x
        ls += x
        n += 1.0
        if i == 0:
            moose.connect( soma, 'axial', dendList[i], 'raxial')
        else:
            moose.connect( dendList[i-1], 'axial', dendList[i], 'raxial')
    return elecid

# def makePassiveSoma( name, length, diameter):
#     numComp = 4
#     x = length/numComp
#     elecid = moose.Neuron( '/library/' + name )
#     soma = moose.vec( elecid.path + '/soma', n=numComp, dtype='Compartment' )

#     ls = x
#     n = 2.0
#     for ind,dend in enumerate(soma):
#         dend.diameter = diameter/n
#         dend.x0 = ls
#         dend.x = ls + x
#         dend.length = x
#         ls += x
#         n += 1.0
#         print( "dend path: ", dend.path )

#     return elecid


    
def extractParms( tab, preStim, endStim, dt ):
    #aoc = sum(tab) - min(tab) * len( tab )
    #aoc = sum(tab) - tab[0] * len( tab )
    start = int( preStim / dt )
    end = int( endStim / dt )
    temp = tab[start:end]
    if displayMoogli:
        aoc = 0.0
        peak = 0.0
    else:
        aoc = sum(temp) - tab[start-1] * len( temp )
        peak = max( temp )
    return [aoc, peak] 

def runTrial( diffusionLength, v, dist, blanks, preStim, postStim ):
    settleTime = params['settleTime']
    vel = moose.vec( '/model/chem/dend/vel' )
    vel.nInit = v * diffusionLength
    #B = moose.vec( '/model/chem/dend/B' )
    #B.nInit = 0.004
    #print 'runTrial(', v, dist, blanks, ')'
    runtime = preStim + postStim + (dist+ blanks*2)/float(v)

    moose.reinit()
    if not displayMoogli:
        moose.start( runtime )
    tabs = moose.vec( '/model/graphs/plot0' )
    #print len(tabs), len( tabs[0].vector )
    stride = int(dist) / numSpine
    #print 'tab number = ', blanks, blanks + dist/2, blanks + dist - 1
    ret = []
    ps = preStim - 4 * params[ 'stimWidth' ]
    for i in range( numSpine ):
        data = tabs[(blanks + i*stride)].vector
        ret.extend( extractParms( data, ps, runtime + settleTime - postStim, tabs[0].dt ) )
    #print ret
    aoc = 0.0
    peak = 0.0
    spatialSum = np.array( [0.0, 0.0] ) # aoc, peak
    settleTime = params['settleTime']
    
    for i in tabs:
        data = i.vector
        spatialSum += np.array( extractParms( data, ps, runtime + settleTime - postStim, tabs[0].dt ) )

    ret.extend( [spatialSum[0], spatialSum[1] ] )
    return ret

def writeXML( seqDtRange, drange, labels, allSimOutput ):
    seqList, seqScore = makeSequence( numSpine )
    sh = allSimOutput.shape
    print "SHAPES = ", sh
    print "seqList = ", seqList
    print "seqScore = ", seqScore
    print 'LEN Labels = ', len(labels)

    assert( len(sh) == 4 and 
            sh[0] == len(seqDtRange) and 
            sh[1] == len(drange) and 
            sh[2] == len(seqList) and 
            sh[3] == len(labels) )
    root = ET.Element( 'twoDplot' )
    yaxis = ET.SubElement( root, 'yaxis' )
    yaxis.set( "title", "seqDt" )
    yaxis.text = ''.join( str(seqDt) + ' ' for seqDt in seqDtRange ) + '\n'
    xaxis = ET.SubElement( root, 'xaxis' )
    xaxis.set( "title", "distance" )
    xaxis.text = ''.join( str(d) + ' ' for d in drange ) + '\n'
    parameters = ET.SubElement( root, 'parameters' )
    p = []
    for j in params:
        p.append( ET.SubElement( parameters, j ) )
        p[-1].text = str( params[j] )

    seqOutput = ET.SubElement( root, 'seqOutput' )
    for iseq in range( len( seqList ) ):
        seqData = ET.SubElement( seqOutput, 'stimSequence' )
        seqData.set( 'index', str( iseq ) )
        seqData.set( 'score', str( seqScore[iseq] ) )
        seqData.text = ''.join( str(j) + ' ' for j in seqList[iseq] )

    simOutput = ET.SubElement( root, 'simOutput' )
    for idt in range( len( seqDtRange ) ):
        dtData = ET.SubElement( simOutput, 'dtData' )
        dtData.set( 'dt', str( seqDtRange[idt]) )
        for idistance in range( len( drange ) ):
            distanceData = ET.SubElement( dtData, 'distanceData' )
            distanceData.set( 'distance', str( drange[idistance] ) )
            for ilab in range( len( labels ) ):
                labelData = ET.SubElement( distanceData, 'labelData' )
                labelData.set( 'title', labels[ilab] )
                y = allSimOutput[idt,idistance,:,ilab]
                labelData.text = ''.join( str(j) + ' ' for j in y )

    tree = ET.ElementTree( root )
    # tree.write( outdir + fname.replace('.p','') + '.' + str( int(params['fnumber']) ) +  '.xml' )
    tree.write( fname.replace('.p','') +  '.xml' )

def convertSeq( arg ):
    x = int( arg )
    ret = [0] * numSpine
    for i in range( numSpine ):
        ret[-i - 1 ] = x % 10
        x /= 10
    return [ret,]

def makeSequence( numSpines ):
    x = range( numSpines )
    a = list( itertools.permutations( x ) )
    skipWeight = 0.64

    linearity = []
    allSlopes = []
    b = []
    #patterns = [ [1,2,3,4,0],[1,4,2,0,3],[2,0,4,3,1],[4,3,2,1,0],[0,1,2,3,4]]
    #patterns = a[0::2]
    if displayMoogli:
        patterns = convertSeq( moogliPattern )
    else:
        patterns = a[0::1]
        # patterns = a[0::10080]
        #patterns = [[0,1,2,3,4],[4,0,3,2,1],[4,3,2,1,0]]

    for i in patterns:   # Do every other sample.
        slope, intercept, r_value, p_value, std_err = stats.linregress(x,i)
        '''
        if slope > 0.1:
            b.append( i )
            linearity.append( r_value**2 )
            #print len( b ), i, linearity[-1]
        '''
        b.append( i )
        #linearity.append( np.sign( slope ) * r_value**2 )
        linearity.append( slope * r_value**2 )
    return b, linearity

def scanDistDtGrid():
    numCompts = 22
    comptLength = 1
    startPos = numCompts * comptLength 
    seqList, seqScore = makeSequence( numSpine )
    #moose.setClock(18, 0.02 )
    #print "DT = ", moose.element( '/model/graphs/plot0').dt
    seqDtRange = (1.0, 2.0, 3.0, 4.0 ) # Interval between stimuli, in sec
    drange = ( 5, 10, 15, 20) # Distance covered, in diffLen.
    # seqDtRange = ( 1.0, 2.0, 3.0)
    # drange =  ( 15, 20, )
    if displayMoogli:
        seqDtRange = ( moogliDt, )
        drange = ( moogliDistance, )
    allSimOutput = []

    ampl = moose.vec( '/model/chem/dend/ampl' )
    print("ampl: ", len(ampl))
    ampl.nInit = params['stimAmplitude']
    blanks = params['blankVoxelsAtEnd']
    preStim = params['preStimTime']
    postStim = params['postStimTime']
    diffusionLength = params['diffusionLength']

    f1 = open('./stimPnts.txt', 'a')
    
    for seqDt in seqDtRange:
        temp2 = []
        temp3 = []
        print >>f1, "dt: ", seqDt
        for d in drange:
            print >>f1, "\td: ", d
            # compt = moose.Compartment( '/model/chem/soma')
            # compt.diameter = params['dendDiameter']/8 
            dend = moose.element( '/model/chem/dend')
            # print("dend diameter: ", dend.diameter)
            moose.le( '/model/chem/dend' )

            phase = moose.vec( '/model/chem/dend/phase' )
            ampl = moose.vec( '/model/chem/dend/ampl' )

            # Old version. Assumes we activate density * d compartments.
            #ampl.nInit = float(v)/float(d)  

            ampl.nInit = params['stimAmplitude']
            Z = moose.vec( '/model/chem/dend/Z' )
            print("Z length: ", len(Z))
            stride = int(d) / numSpine
            Z.nInit = 0
            phase.nInit = 10000.0

            temp = []
            slopes = []
            pos = range(startPos+1,len(Z),2)
            print pos
            for seq in seqList:
                print '.',
                sys.stdout.flush()
                for j in range( numSpine ):
                    #k = (blanks + j * stride)
                    k = pos[j*stride]
                    print 'k is', k
                    print >>f1, '\t\tk: ', k 
                    Z[ k ].nInit = 1
                    print 'volumes are',Z[k-2].volume,Z[k-1].volume,Z[k].volume,Z[k+1].volume,Z[k+2].volume
                    phase[ k ].nInit = preStim + seq[j] * seqDt
                temp.append( runTrial( diffusionLength, seqDt, d, blanks, preStim, postStim ))
                #print temp
            print >>f1, "\n\n"
            simOut = np.array( temp )
            temp3.append( simOut )
            print seqDt, d #temp[-1], temp[-2], 
        allSimOutput.append( temp3 )
    return allSimOutput, seqDtRange, drange


def main():
    global displayMoogli
    global moogliPattern
    global moogliDt
    global moogliDistance
    # global rdes

    # if len( sys.argv ) == 2:
    #     print "Usage: ", sys.argv[0], " --param value --moogli distance dt sequence"
    #     quit()

    displayMoogli = False
    for ii in range( len( sys.argv ) ):
        if sys.argv[ii][:2] == '--':
            argName = sys.argv[ii][2:]
            if argName in params:
                params[argName] *= float(sys.argv[ii+1])
            if argName == 'moogli':
                displayMoogli = True
                moogliDistance = float( sys.argv[ii+1] )
                moogliDt = float( sys.argv[ii+2] )
                moogliPattern = sys.argv[ii+3]

    # moose.le('/library')
    # moose.le('/library/soma')
                
    # moose.le( '/model/chem' )

    diffusionLength = params['diffusionLength']
    dendLength = params['dendLength']    #100 microns
    diffusionLength = params['diffusionLength']   #1 micron
    library = moose.Neutral( '/library' )
    fac = 1.0
#def makeChemProto( name, stimAmpl = 1, diffLength = 1e-6, ):
    makeChemProto( 'spread', stimAmpl = 1,
            diffLength = diffusionLength )
    makePassiveSoma( 'cell', params['dendLength'] * fac ,  params['dendDiameter'] )

    rdes = rd.rdesigneur(
        useGssa = False,
        turnOffElec = True,
        chemPlotDt = 0.02,
        #chemDt = 0.01,
        #diffDt = 0.01,
        diffusionLength = diffusionLength,
        cellProto = [['./flex_tap.p', 'elec']],
	    #cellProto = [['./models/' + fname, 'soma']],	# to take filenames as arguments
        chemProto = [['spread', 'spread']],
        chemDistrib = [['spread', '#', 'install', '1' ]],
        plotList = [['#', '1', 'dend/A', 'n', '# of A'],
            ['#', '1', 'dend/B', 'n', '# of B'],
            ['#', '1', 'dend/Ca', 'n', '# of Ca']
            ],
        # moogList = [
        #     ['soma', '1', 'dend/A', 'n', 'A n', 0, 4],
        #     ['soma', '1', 'dend/B', 'n', 'B n', 0, 2.5],
            # ['soma', '1', 'dend/Ca', 'n', 'Ca n', 0.0, 1]
        # ]
         moogList = [['#', 'x>5e-6 && x<100e-6', '.', 'Vm', 'Vm'],
            # ['#', 'y<1e-6', 'dend/A', 'n', 'A n', 0, 4],
            # ['#', 'y<1e-6', 'dend/B', 'n', 'B n', 0, 2.5],
            # ['#', 'x>5e-6 && x<100e-6', 'dend/Ca', 'n', 'Ca n', 0.0, 1]

        ]
    )
    rdes.buildModel()

    # return

    allSimOutput, seqDtRange, drange = scanDistDtGrid()

    if displayMoogli:
        preStim = params['preStimTime']
        postStim = params['postStimTime']
        blanks = params['blankVoxelsAtEnd']
        runtime =  postStim + preStim + (drange[0]+ blanks*2)/float(seqDtRange[0])
        rdes.displayMoogli( 0.1, runtime, 0.0 )
    else:
        labels = ['aoc0', 'peak0', 'aoc1', 'peak1', 'aoc2', 'peak2', 'aoc3', 'peak3', 'aoc4', 'peak4', 'aocSum', 'peakSum' ]
        # labels = ['aoc0', 'peak0', 'aoc1', 'peak1', 'aoc2', 'peak2', 'aoc3', 'peak3', 'aoc4', 'peak4', 'aoc5', 'peak5', 'aoc6', 'peak6', 'aoc7', 'peak7',  'aocSum', 'peakSum' ]
        writeXML( seqDtRange, drange, labels, np.array( allSimOutput) )
    
    # with open('fileList.txt', 'a') as f:
        # f.write("\nDetails: \nDistance: %s Dt: %s Seq: %s DendLength: %s Blanks: %s " %(moogliDistance, moogliDt, moogliPattern, dendLength*fac, blanks) )

if __name__ == '__main__':
    main()

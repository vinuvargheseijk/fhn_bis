import moose
import rdesigneur as rd


numDendSegments=100
comptLenBuff=14.4e-06
comptLen=0.216e-06
#comptDia=1e-06
RM=1.0
RA=10.0
CM=0.001

def makeChemProto(name='hydra'):
        chemCompt=moose.Neutral('/library/'+name)
        compt=moose.CubeMesh('/library/'+name + '/' + name)
        A=moose.Pool(compt.path+'/A')
        B=moose.Pool(compt.path+'/B')
        Adot = moose.Function( A.path + '/Adot' )
        Bdot = moose.Function( B.path + '/Bdot' )
        #Adot.expr="1*(-x0+x1*(0.05*+(x0^2/(1^2+x0^2))))"
        #Bdot.expr="-1*(-x0+x1*(0.05*+(x0^2/(1^2+x0^2))))"
        Adot.expr="0*x0+0*x1"
        Bdot.expr="0*x0+0*x1"
 
        
        print "$$$$> ", Adot, Bdot
        print Adot.expr, Bdot.expr
        print moose.showmsg(Adot)
        Adot.x.num = 2
        Bdot.x.num = 2
        A.nInit=1
        B.nInit=1
        moose.connect( A, 'nOut', Adot.x[0], 'input' )
        moose.connect( B, 'nOut', Adot.x[1], 'input' )
        moose.connect( Adot, 'valueOut', A, 'increment' )

        moose.connect( A, 'nOut', Bdot.x[0], 'input' )
        moose.connect( B, 'nOut', Bdot.x[1], 'input' )
        moose.connect( Bdot, 'valueOut', B, 'increment' )
        return compt

def makeDendProto():
    dend=moose.Neuron('/library/dend')
    prev=rd.buildCompt(dend,'soma',RM=RM,RA=RA,CM=CM,dia=0.3e-06,x=0,dx=comptLenBuff)
    x=comptLenBuff
    y=0.0
    comptDia=0.3e-06

    for i in range(numDendSegments):
      dx=comptLen
      dy=0
      #comptDia +=1.7e-08
      compt=rd.buildCompt(dend,'dend'+str(i),RM=RM,RA=RA,CM=CM,x=x,y=y,dx=dx,dy=dy,dia=comptDia)
      moose.connect(prev,'axial',compt,'raxial')
      prev=compt
      x+=dx
      y+=dy
      
    compt=rd.buildCompt(dend,'dendL',RM=RM,RA=RA,CM=CM,x=x,y=y,dx=comptLenBuff,dy=dy,dia=comptDia)
    moose.connect(prev,'axial',compt,'raxial')
    return dend


library=moose.Neutral('/library')        
makeChemProto()
makeDendProto()
rdes = rd.rdesigneur(
    cellProto = [['./flex_tap.p', 'elec']],
    chemProto = [['hydra','chem']],
    chemDistrib = [['chem','#','install','1']],
    #stimList = [['#spine#', '1', '.', 'inject', 't * 25e-9' ]],
    stimList = [['#', '1', '.', 'inject', 't * 25e-9' ]],
    #plotList = [['#', '1', '.', 'Vm', 'Membrane potential'],
     #   ['#', '1', 'Ca_conc', 'Ca', 'Ca conc (uM)']],
    moogList = [['#', '1', '.', 'Vm', 'Soma potential']]
)

rdes.buildModel()

moose.reinit()
rdes.displayMoogli( 0.0002, 0.1 )
